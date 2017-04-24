#include "defs.hpp"

#include "node_server.hpp"
#include "node_client.hpp"
#include "future.hpp"
#include "options.hpp"

#include <chrono>
#include <string>
#include <utility>
#include <vector>

#include <fenv.h>
#if !defined(_MSC_VER)
#include <unistd.h>
#else
#include <float.h>
#endif

#include <hpx/hpx_init.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/lcos/broadcast.hpp>

options opts;

void initialize(options _opts, std::vector<hpx::id_type> const& localities) {
	options::all_localities = localities;
	opts = _opts;
#if !defined(_MSC_VER)
	feenableexcept(FE_DIVBYZERO);
	feenableexcept(FE_INVALID);
	feenableexcept(FE_OVERFLOW);
#else
	_controlfp(_EM_INEXACT | _EM_DENORMAL | _EM_INVALID, _MCW_EM);
#endif


}

HPX_PLAIN_ACTION(initialize, initialize_action);
HPX_REGISTER_BROADCAST_ACTION_DECLARATION (initialize_action)
HPX_REGISTER_BROADCAST_ACTION (initialize_action)

int hpx_main(int argc, char* argv[]) {
	printf("###########################################################\n");
#if defined(__AVX512F__)
	printf("Compiled for AVX512 SIMD architectures.\n");
#elif defined(__AVX2__)
	printf("Compiled for AVX2 SIMD architectures.\n");
#elif defined(__AVX__)
	printf("Compiled for AVX SIMD architectures.\n");
#elif defined(__SSE2__ )
	printf("Compiled for SSE2 SIMD architectures.\n");
#else
	printf("Not compiled for a known SIMD architecture.\n");
#endif
	printf("###########################################################\n");

	printf("Running\n");

	try {
		if (opts.process_options(argc, argv)) {
			auto all_locs = hpx::find_all_localities();
			hpx::lcos::broadcast < initialize_action
					> (all_locs, opts, all_locs).get();

			node_client root_id = hpx::new_ < node_server > (hpx::find_here());
			node_client root_client(root_id);
			node_server* root = root_client.get_ptr().get();

			int ngrids = 0;
			if (opts.found_restart_file) {

				//	set_problem(null_problem);

				const std::string fname = opts.restart_filename;
				printf("Loading from %s...\n", fname.c_str());
				if (opts.output_only) {
					const std::string oname = opts.output_filename;
					root->load_from_file_and_output(fname, oname,
							opts.data_dir);
				} else {
					root->load_from_file(fname, opts.data_dir);
					ngrids = root->regrid(root_client.get_gid(), true);

				}
				printf("Done. \n");
			} else {
				for (integer l = 0; l < opts.max_level; ++l) {
					ngrids = root->regrid(root_client.get_gid(), false);
					printf("---------------Created Level %i---------------\n\n",
							int(l + 1));
				}
				ngrids = root->regrid(root_client.get_gid(), false);
				printf("---------------Regridded Level %i---------------\n\n",
						int(opts.max_level));
			}

			if (!opts.output_only) {
				hpx::async(&node_server::amr_driver, root).get();
			}
		}
	} catch (...) {
		throw;
	}
	printf("Exiting...\n");
	return hpx::finalize();
}

int main(int argc, char* argv[]) {
	std::vector<std::string> cfg = { "hpx.commandline.allow_unknown=1", // HPX should not complain about unknown command line options
			"hpx.scheduler=local-priority-lifo", // use LIFO scheduler by default
			"hpx.parcel.mpi.zero_copy_optimization!=0" // Disable the usage of zero copy optimization for MPI...
			};

	hpx::register_pre_shutdown_function(
			[]() {options::all_localities.clear();});

	hpx::init(argc, argv, cfg);
}
