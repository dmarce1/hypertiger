/*
 * options.hpp
 *
 *  Created on: Nov 13, 2015
 *      Author: dmarce1
 */

#ifndef OPTIONS_HPP_
#define OPTIONS_HPP_

#include <string>
#include <vector>
#include <hpx/hpx.hpp>
#include "defs.hpp"


class options {

	std::string exe_name;

	bool cmp(const char* str1, const char* str2);
	bool cmp(const std::string str1, const char* str2);
	void show_help();
public:
	integer max_level;
	real xscale;
	std::string restart_filename;
	bool found_restart_file;
	std::string output_filename;
	bool output_only;
	real output_dt;
	real stop_time;
    integer stop_step;
    bool disable_output;
    bool parallel_silo;
    bool silo_planes_only;
    std::string data_dir;

	template<class Arc>
	void serialize(Arc& arc, unsigned) {
		arc & parallel_silo;
		arc & silo_planes_only;
		arc & stop_time;
		arc & max_level;
		arc & xscale;
		arc & restart_filename;
		arc & found_restart_file;
		arc & output_filename;
		arc & output_only;
		arc & output_dt;
        arc & stop_step;
        arc & disable_output;
        arc & data_dir;
	}

	bool process_options(int argc, char* argv[]);

    static std::vector<hpx::id_type> all_localities;
};

#ifndef IN_OPTIONS_CPP
extern options opts;
#endif

#endif /* OPTIONS_HPP_ */
