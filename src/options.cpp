/*
 * options.cpp
 *
 *  Created on: Nov 13, 2015
 *      Author: dmarce1
 */

#include "defs.hpp"
#define IN_OPTIONS_CPP
#include "options.hpp"
#include <math.h>
#include "grid.hpp"

#define HELP_OPT "-Help"
#define PROBLEM_OPT "-Problem"
#define RESTART_OPT "-Restart"
#define OUTPUT_OPT "-Output"
#define DATA_DIR_OPT "-Datadir"
#define XSCALE_OPT "-Xscale"
#define ODT_OPT "-Odt"
#define DISABLEOUTPUT_OPT "-Disableoutput"
#define STOPTIME_OPT "-Stoptime"
#define STOPSTEP_OPT "-Stopstep"
#define PARALLEL_SILO_OPT "-ParallelSilo"
#define SILO_PLANES_ONLY_OPT "-SiloPlanesOnly"

#define MAX_LEVEL_OPT "-Max_level"

bool options::cmp(const char* str1, const char* str2) {
	return strncmp(str1, str2, strlen(str2)) == 0;
}

bool options::cmp(const std::string str1, const char* str2) {
	return strncmp(str1.c_str(), str2, strlen(str2)) == 0;
}

void options::show_help() {
}

bool options::process_options(int argc, char* argv[]) {
	bool rc;
	rc = true;
	parallel_silo = false;
	silo_planes_only = false;
	max_level = 3;
	found_restart_file = false;
	output_only = false;
	xscale = 1.0;
	exe_name = std::string(argv[0]);
	output_dt = -1;
	stop_time = std::numeric_limits<real>::max() - 1;
	stop_step = std::numeric_limits<integer>::max() / 10;
	disable_output = false;
	bool vomega_found = false;

	for (integer i = 1; i < argc; ++i) {
		if (cmp(argv[i], HELP_OPT)) {
			rc = false;
		} else if (cmp(argv[i], PROBLEM_OPT)) {
			std::string prob(argv[i] + strlen(PROBLEM_OPT) + 1);
		} else if (cmp(argv[i], RESTART_OPT)) {
			restart_filename = std::string(argv[i] + strlen(RESTART_OPT) + 1);
			found_restart_file = true;
		} else if (cmp(argv[i], DATA_DIR_OPT)) {
			data_dir = std::string(argv[i] + strlen(DATA_DIR_OPT) + 1);
			data_dir += "/";
		} else if (cmp(argv[i], PARALLEL_SILO_OPT)) {
			parallel_silo = true;
		} else if (cmp(argv[i], SILO_PLANES_ONLY_OPT)) {
			silo_planes_only = true;
		} else if (cmp(argv[i], OUTPUT_OPT)) {
			output_filename = std::string(argv[i] + strlen(OUTPUT_OPT) + 1);
			output_only = true;
		} else if (cmp(argv[i], MAX_LEVEL_OPT)) {
			max_level = atoi(argv[i] + strlen(MAX_LEVEL_OPT) + 1);
		} else if (cmp(argv[i], XSCALE_OPT)) {
			xscale = atof(argv[i] + strlen(XSCALE_OPT) + 1);
		} else if (cmp(argv[i], ODT_OPT)) {
			output_dt = atof(argv[i] + strlen(ODT_OPT) + 1);
		} else if (cmp(argv[i], DISABLEOUTPUT_OPT)) {
			disable_output = true;
		} else if (cmp(argv[i], STOPTIME_OPT)) {
			stop_time = atof(argv[i] + strlen(STOPTIME_OPT) + 1);
		} else if (cmp(argv[i], STOPSTEP_OPT)) {
			stop_step = atoi(argv[i] + strlen(STOPSTEP_OPT) + 1);
		} else {
			printf("Unknown option - %s\n", argv[i]);
			abort();
		}
	}
	if (output_dt < 0) {
		output_dt = 1.0e-2;
	}
	if (!rc) {
		show_help();
	}
	return rc;
}

std::vector<hpx::id_type> options::all_localities = { };
