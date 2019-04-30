/*****************************************************************************\
 *  site_factor_xfactor.c - xfactor site_factor plugin
*****************************************************************************
 *  Copyright (C) 2019 SchedMD LLC
 *  Written by Alejandro Sanchez <alex@schedmd.com>
 *
 *  This file is part of Slurm, a resource management program.
 *  For details, see <https://slurm.schedmd.com/>.
 *  Please also read the included file: DISCLAIMER.
 *
 *  Slurm is free software; you can redistribute it and/or modify it under
 *  the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  In addition, as a special exception, the copyright holders give permission
 *  to link the code of portions of this program with the OpenSSL library under
 *  certain conditions as described in each individual source file, and
 *  distribute linked combinations including the two. You must obey the GNU
 *  General Public License in all respects for all of the code used other than
 *  OpenSSL. If you modify file(s) with this exception, you may extend this
 *  exception to your version of the file(s), but you are not obligated to do
 *  so. If you do not wish to do so, delete this exception statement from your
 *  version.  If you delete this exception statement from all source files in
 *  the program, then also delete it here.
 *
 *  Slurm is distributed in the hope that it will be useful, but WITHOUT ANY
 *  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 *  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 *  details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with Slurm; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301  USA.
\*****************************************************************************/

#define _GNU_SOURCE

#include "src/common/slurm_xlator.h"
#include "src/common/log.h"

#include <math.h>

/*
 * These variables are required by the generic plugin interface.  If they
 * are not found in the plugin, the plugin loader will ignore it.
 *
 * plugin_name - A string giving a human-readable description of the
 * plugin.  There is no maximum length, but the symbol must refer to
 * a valid string.
 *
 * plugin_type - A string suggesting the type of the plugin or its
 * applicability to a particular form of data or method of data handling.
 * If the low-level plugin API is used, the contents of this string are
 * unimportant and may be anything.  Slurm uses the higher-level plugin
 * interface which requires this string to be of the form
 *
 *	<application>/<method>
 *
 * where <application> is a description of the intended application of
 * the plugin (e.g., "auth" for Slurm authentication) and <method> is a
 * description of how this plugin satisfies that application.  Slurm will
 * only load authentication plugins if the plugin_type string has a prefix
 * of "auth/".
 *
 * plugin_version - an unsigned 32-bit integer containing the Slurm version
 * (major.minor.micro combined into a single number).
 */
const char	*plugin_name		= "xfactor site_factor plugin";
const char	*plugin_type		= "site_factor/xfactor";
const uint32_t	plugin_version		= SLURM_VERSION_NUMBER;

uint32_t xfactor_min_time = 1;		/* Minimum time limit in minutes. */
uint32_t xfactor_max = NICE_OFFSET;	/* Maximum weightened xfactor value. */
uint32_t xfactor_weight = 1;		/* Weight for xfactor factor. */

static void _parse_parameters(void)
{
	char *params = NULL, *tmp_ptr = NULL;
	uint32_t weight = 0, min_time = 0, max = 0;

	params = slurm_get_priority_site_factor_params();

	if (!params) {
		error("%s: PrioritySiteFactorParameters not set.", plugin_type);
		return;
	}

	tmp_ptr = xstrcasestr(params, "xfactor_min_time=");
	if (!tmp_ptr) {
		error("%s: xfactor_min_time not configured.", plugin_type);
		goto cleanup;
	}
	min_time = atoi(tmp_ptr + 17);
	if (min_time < 1 || min_time > 129600) {
		error("%s: invalid xfactor_min_time value.", plugin_type);
		goto cleanup;
	}
	xfactor_min_time = min_time;

	tmp_ptr = xstrcasestr(params, "xfactor_max=");
	if (!tmp_ptr) {
		error("%s: xfactor_max not configured.", plugin_type);
		goto cleanup;
	}
	max = atoi(tmp_ptr + 12);
	if (max < 1 || max > NICE_OFFSET) {
		error("%s: invalid xfactor_max value.", plugin_type);
		goto cleanup;
	}
	xfactor_max = max;

	tmp_ptr = xstrcasestr(params, "xfactor_weight=");
	if (!tmp_ptr) {
		error("%s: xfactor_weight not configured.", plugin_type);
		goto cleanup;
	}
	weight = atoi(tmp_ptr + 15);
	if (weight < 1 || weight > NICE_OFFSET) {
		error("%s: invalid xfactor_weight value.", plugin_type);
		goto cleanup;
	}
	xfactor_weight = weight;

cleanup:
	xfree(params);
	if (slurm_get_debug_flags() & DEBUG_FLAG_PRIO)
		debug("%s: xfactor_min_time=%u, xfactor_max=%u, xfactor_weight=%u",
		      plugin_type, xfactor_min_time, xfactor_max,
		      xfactor_weight);

	return;
}

extern int init(void)
{
	debug("%s: %s loaded", __func__, plugin_name);
	_parse_parameters();

	return SLURM_SUCCESS;
}

extern int fini(void)
{
	debug("%s: unloading %s", __func__, plugin_name);

	return SLURM_SUCCESS;
}

extern void site_factor_p_reconfig(void)
{
	 _parse_parameters();

	return;
}

static uint32_t _calc_factor(struct job_record *job_ptr)
{
	uint32_t factor = 0, delta = 0, tmp_time = 1;
	uint64_t weightened_factor = 0;
	double quotient = 0.0;
	time_t now = time(NULL);

	if (!xfactor_weight || !job_ptr->details ||
	    !job_ptr->details->accrue_time)
		return factor;

	if (now > job_ptr->details->accrue_time)
		delta = now - job_ptr->details->accrue_time;

	if (!delta)
		return factor;

	if (job_ptr->time_limit != NO_VAL)
		tmp_time = job_ptr->time_limit;
	else if (job_ptr->part_ptr &&
		 (job_ptr->part_ptr->max_time != INFINITE))
		tmp_time = job_ptr->part_ptr->max_time;

	tmp_time = MAX(tmp_time, xfactor_min_time);
	quotient = (double)delta / (double)tmp_time;
	factor = lround(quotient);

	weightened_factor = factor * xfactor_weight;
	factor = MIN(weightened_factor, xfactor_max);

	if (slurm_get_debug_flags() & DEBUG_FLAG_PRIO)
		debug("%s: weightened site_factor=%u", plugin_type, factor);

	return factor;
}

extern void site_factor_p_set(struct job_record *job_ptr)
{
	job_ptr->site_factor = _calc_factor(job_ptr) + NICE_OFFSET;

	return;
}

static int _update(void *x, void *ignored)
{
	struct job_record *job_ptr = (struct job_record *) x;

	if (IS_JOB_PENDING(job_ptr))
		job_ptr->site_factor = _calc_factor(job_ptr) + NICE_OFFSET;

	return SLURM_SUCCESS;
}

extern void site_factor_p_update(void)
{
	list_for_each(job_list, _update, NULL);

	return;
}
