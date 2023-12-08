/*
 *
 * Copyright 2021-2023 Software Radio Systems Limited
 *
 * This file is part of srsRAN.
 *
 * srsRAN is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * srsRAN is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * A copy of the GNU Affero General Public License can be found in
 * the LICENSE file in the top-level directory of this distribution
 * and at http://www.gnu.org/licenses/.
 *
 */

#include "srsran/phy/upper/channel_coding/channel_coding_factories.h"
#include "srsran/support/benchmark_utils.h"
#include "srsran/support/srsran_test.h"
#include <getopt.h>
#include <random>
#include <iostream>
#include <cstdint>
#include <climits>
#include <chrono>
#include <vector>
#include <algorithm>
#include <fstream> 
#include <iomanip>
//#include </home/ubuntu/perf_cycle_codes/pmuv3.h>
//#include "/home/ubuntu/perf_cycle_codes/c_pmuv3_functions.h"

/* HEADERS FOR PMUV3 */
// SPDX-License-Identifier: GPL-2.0

extern "C"{

#include "/home/ubuntu/perf_cycle_codes/pmuv3_plugin.h"
//#include "/home/ubuntu/perf_cycle_codes/pmuv3.h"

}

extern "C"{

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <linux/perf_event.h>
#include <linux/kernel.h>
//#include <perf/cpumap.h>
#include </home/ubuntu/linux-5.17.4/tools/lib/perf/include/perf/cpumap.h>
#include </home/ubuntu/linux-5.17.4/tools/lib/perf/include/perf/threadmap.h> //Modified for LDPC
//#include <threadmap.h>
//#include <perf/threadmap.h>
//#include <evsel.h>
//#include <perf/evsel.h>
//#include <internal/evsel.h>
//#include <evsel.h>
#include </home/ubuntu/linux-5.17.4/tools/lib/perf/include/perf/evsel.h>
//#include <internal/tests.h>
#include </home/ubuntu/linux-5.17.4/tools/lib/perf/include/internal/tests.h>
//#include "tests.h"
#include <sys/mman.h>
#include <time.h>
#include <inttypes.h> // Include inttypes.h for PRIu64 macro
#include <asm/unistd.h> 
}
/*25 MHZ timer generic */
#include <stdio.h>
#include <stdint.h>
#include <time.h> // For time measurements in C
#include <unistd.h> // For sleep functionality in C 
/********************************************************** HEADERS FOR PMUV3 ************************************************************/
//#define __printf(a, b)  __attribute__((format(printf, a, b)))

#if 0
//******************HEADERS FOR LIBPERF DEFINITIONS*****************************************************
//#include <linux/compiler.h>
#include </home/ubuntu/linux-5.17.4/tools/include/linux/compiler.h>
//#include <perf/core.h>
#include </home/ubuntu/linux-5.17.4/tools/lib/perf/include/perf/core.h>
#include </home/ubuntu/linux-5.17.4/tools/lib/perf/include/internal/lib.h>
//#include <internal/lib.h>
//#include "internal.h"
#include </home/ubuntu/linux-5.17.4/tools/lib/perf/internal.h>
#endif
//*******************************HEADERS FOR NEW DUMMY THREADMAP.C########################
//#include <asm/bug.h>
//#include </home/ubuntu/linux-5.17.4/tools/include/asm/bug.h>
//#include </home/ubuntu/linux-5.17.4/tools/lib/perf/include/internal/threadmap.h>

/*****************************************************************************************************************************************
static int __base_pr(enum libperf_print_level level __maybe_unused, const char *format,
                     va_list args)
{
        return vfprintf(stderr, format, args);
}

static libperf_print_fn_t __libperf_pr = __base_pr;

//__printf(2, 3)
#if 0
void libperf_print(enum libperf_print_level level, const char *format, ...)
{
        va_list args;

        if (!__libperf_pr)
                return;

        va_start(args, format);
        __libperf_pr(level, format, args);
        va_end(args);
}
#endif





extern "C" {

static int libperf_print(enum libperf_print_level level,
                         const char *fmt, va_list ap)
{
        //return 0;
        return vfprintf(stderr, fmt, ap);
} 
}
void libperf_init(libperf_print_fn_t fn) 
{
        page_size = sysconf(_SC_PAGE_SIZE);
        __libperf_pr = fn; 
}

void perf_thread_map__set_pid(struct perf_thread_map *map, int thread, pid_t pid)
{
        map->map[thread].pid = pid;
}

struct perf_thread_map *perf_thread_map__new_dummy(void)
{
        struct perf_thread_map *threads = thread_map__alloc(1);

        if (threads != NULL) {
                perf_thread_map__set_pid(threads, 0, -1);
                threads->nr = 1;
                refcount_set(&threads->refcnt, 1); 
        }
        return threads;
}
******************************************************************************************************************************/

/*
static int test_stat_user_read(int event)
{
        struct perf_counts_values counts = { .val = 0 };
        struct perf_thread_map *threads;
        struct perf_evsel *evsel;
        struct perf_event_mmap_page *pc;
        struct perf_event_attr attr = {
                .type   = PERF_TYPE_HARDWARE,
                .config = event,

                .config1 = 0x2,         // Request user access 
#ifdef __aarch64__
                .config1 = 0x2,         // Request user access 
#endif 
        };
        int err, i;

        __u64 start, end, last = 0;

        threads = perf_thread_map__new_dummy();
        __T("failed to create threads", threads);

        perf_thread_map__set_pid(threads, 0, 0);

        evsel = perf_evsel__new(&attr);
        __T("failed to create evsel", evsel);

        err = perf_evsel__open(evsel, NULL, threads);
        __T("failed to open evsel", err == 0);

        err = perf_evsel__mmap(evsel, 0);
        __T("failed to mmap evsel", err == 0);

        pc = perf_evsel__mmap_base(evsel, 0, 0);
        __T("failed to get mmapped address", pc);

#if 1
//volatile int iter = 1000000000;
volatile int iter = 1000;
perf_evsel__read(evsel, 0, 0, &counts);
start = counts.val;
while (iter--) {
        p[random() & 0xFFFF] = random(); // Writing to the array with random numbers
    }
perf_evsel__read(evsel, 0, 0, &counts);
end = counts.val;
printf("Cycle count: %llu\n", end - start);
#endif
        perf_evsel__munmap(evsel);
        perf_evsel__close(evsel);
        perf_evsel__delete(evsel);

        perf_thread_map__put(threads);
        return 0;
}

int test_evsel(int argc, char **argv)
{
        __T_START;

        libperf_init(libperf_print);

        //test_stat_cpu();
        //test_stat_thread();
        //test_stat_thread_enable();
        test_stat_user_read(PERF_COUNT_HW_INSTRUCTIONS);
//      test_stat_user_read(PERF_COUNT_SW_CPU_CLOCK);
        test_stat_user_read(PERF_COUNT_HW_CACHE_L1D);
        test_stat_user_read(PERF_COUNT_HW_CPU_CYCLES);
        test_stat_user_read(PERF_COUNT_HW_CACHE_REFERENCES);
        test_stat_user_read(PERF_COUNT_HW_CACHE_MISSES);
//      test_stat_user_read(PERF_COUNT_HW_BRANCH_INSTRUCTIONS);
        test_stat_user_read(PERF_COUNT_HW_BUS_CYCLES);
        test_stat_user_read(PERF_COUNT_HW_BRANCH_MISSES);
        test_stat_user_read(PERF_COUNT_HW_STALLED_CYCLES_FRONTEND);
        test_stat_user_read(PERF_COUNT_HW_STALLED_CYCLES_BACKEND);
//      test_stat_user_read(PERF_COUNT_HW_REF_CPU_CYCLES);
        //test_stat_user_read(PERF_COUNT_HW_REF_CPU_CYCLES);
        //test_stat_read_format();

        __T_END;
        return tests_failed == 0 ? 0 : -1;
}


int tests_failed;
int tests_verbose;
int main(int argc, char **argv)
{
__T("test evsel", !test_evsel(argc, argv));
        return 0;
}

*/

//uint64_t read_timer() {

   // uint64_t value;
    // Read Physical Timer (Cycle) Counter register
    // This requires user space access to be enabled in cntkctl_el1
   // __asm __volatile("dsb sy");
   // __asm __volatile("MRS %0, PMCCNTR_EL0" : "=r"(value) : : "memory");
   // __asm __volatile("dsb sy");
   // return value;
//}
// clock
#define GENERATE_READ_SYSTEM_REGISTER(ResultType, FuncName, Reg)               \
  inline ResultType FuncName() {                                               \
    uint64_t Res;                                                             \
    __asm__ volatile("mrs \t%0," #Reg : "=r"(Res));                            \
    return Res;                                                               \
  }

GENERATE_READ_SYSTEM_REGISTER(uint64_t, readCycleCount,
                              cntvct_el0)

namespace {
std::mt19937 rgen(0);
std::string  enc_type        = "neon";
unsigned     nof_repetitions = 1000;
bool         silent          = false;
} // namespace

static void usage(const char* prog)
{
  fmt::print("Usage: {} [-R repetitions] [-s silent]\n", prog);
  fmt::print("\t-R Repetitions [Default {}]\n", nof_repetitions);
  fmt::print("\t-T Encoder type generic, avx2 or neon [Default {}]\n", enc_type);
  fmt::print("\t-s Toggle silent operation [Default {}]\n", silent);
  fmt::print("\t-h Show this message\n");
}

static void parse_args(int argc, char** argv)
{
  int opt = 0;
  while ((opt = getopt(argc, argv, "R:T:sh")) != -1) {
    switch (opt) {
      case 'R':
        nof_repetitions = std::strtol(optarg, nullptr, 10);
        break;
      case 'T':
        enc_type = std::string(optarg);
        break;
      case 's':
        silent = (!silent);
        break;
      case 'h':
      default:
        usage(argv[0]);
        exit(0);
    }
  }
}


using namespace srsran;
using namespace srsran::ldpc;
double calculatePercentile50(const std::vector<int>& data) {
    if (data.empty()) {
        // Handle the case where the data vector is empty
        return 0.0; // or some default value
    }

    // Make a copy of the data to avoid modifying the original vector
    std::vector<int> sortedData = data;
    std::sort(sortedData.begin(), sortedData.end());

    // Calculate the index corresponding to the 50th percentile
    double index = 0.5 * (sortedData.size() - 1);

    // Calculate the fractional and integer parts of the index
    double fractionalPart = std::fmod(index, 1.0);
    int integerPart = static_cast<int>(index);

    // Interpolate between the two nearest values
    int lowerValue = sortedData[integerPart];
    int upperValue = sortedData[integerPart + 1];
    double interpolatedValue = lowerValue + fractionalPart * (upperValue - lowerValue);

    return interpolatedValue;
}
double calculatePercentile(const std::vector<int>& data, double percentile) {
    if (data.empty()) {
        // Handle the case where the data vector is empty
        return 0.0; // or some default value
    }

    // Make a copy of the data to avoid modifying the original vector
    std::vector<int> sortedData = data;
    std::sort(sortedData.begin(), sortedData.end());

    // Calculate the index corresponding to the desired percentile
    double index = percentile * (sortedData.size() - 1);

    // Calculate the fractional and integer parts of the index
    double fractionalPart = std::fmod(index, 1.0);
    int integerPart = static_cast<int>(index);

    // Interpolate between the two nearest values
    int lowerValue = sortedData[integerPart];
    int upperValue = sortedData[integerPart + 1];
    double interpolatedValue = lowerValue + fractionalPart * (upperValue - lowerValue);

    return interpolatedValue;
}
#if 0
extern "C"{
volatile int rand_arr[0xFFFF]; // Array declaration

void nano_sleep(long nanoseconds) {
    struct timespec req, rem;
    req.tv_sec = 0;
    req.tv_nsec = nanoseconds;

    // Perform the nanosleep
    while (nanosleep(&req, &rem) == -1) {
        req = rem;
    }
}
static int libperf_print(enum libperf_print_level level,
			 const char *fmt, va_list ap)
{
	//return 0;
	return vfprintf(stderr, fmt, ap);
}
}

#endif
#if 0
static int test_stat_user_read(int event)
{
	struct perf_counts_values counts = { .val = 0 };
	struct perf_thread_map *threads;
	struct perf_evsel *evsel;
	struct perf_event_mmap_page *pc;
	struct perf_event_attr attr = {
		.type	= PERF_TYPE_HARDWARE,
		.config	= static_cast<__u64>(event),
		
		.config1 = 0x2,		// Request user access 
	};
	int err;

	__u64 start, end;

	threads = perf_thread_map__new_dummy();
	__T("failed to create threads", threads);

	perf_thread_map__set_pid(threads, 0, 0);

	evsel = perf_evsel__new(&attr);
	__T("failed to create evsel", evsel);

	err = perf_evsel__open(evsel, NULL, threads);
	__T("failed to open evsel", err == 0);

	err = perf_evsel__mmap(evsel, 0);
	__T("failed to mmap evsel", err == 0);

//	pc = perf_evsel__mmap_base(evsel, 0, 0);
	pc = static_cast<perf_event_mmap_page*>(perf_evsel__mmap_base(evsel, 0, 0));

	__T("failed to get mmapped address", pc);



//volatile int iter = 1000000000;
volatile int iter = 1000;
perf_evsel__read(evsel, 0, 0, &counts);
start = counts.val;
while (iter--) {
        rand_arr[random() & 0xFFFF] = random(); // Writing to the array with random numbers
    }
perf_evsel__read(evsel, 0, 0, &counts);
end = counts.val;
printf("Cycle count: %llu\n", end - start);
	perf_evsel__munmap(evsel);
	perf_evsel__close(evsel);
	perf_evsel__delete(evsel);

	perf_thread_map__put(threads);
	return 0;
}

int test_evsel(int argc, char **argv)
{
        __T_START;

        libperf_init(libperf_print);

        //test_stat_cpu();
        //test_stat_thread();
        //test_stat_thread_enable();
        test_stat_user_read(PERF_COUNT_HW_INSTRUCTIONS);
//      test_stat_user_read(PERF_COUNT_SW_CPU_CLOCK);
        test_stat_user_read(PERF_COUNT_HW_CACHE_L1D);
        test_stat_user_read(PERF_COUNT_HW_CPU_CYCLES);
        test_stat_user_read(PERF_COUNT_HW_CACHE_REFERENCES);
        test_stat_user_read(PERF_COUNT_HW_CACHE_MISSES);
//      test_stat_user_read(PERF_COUNT_HW_BRANCH_INSTRUCTIONS);
        test_stat_user_read(PERF_COUNT_HW_BUS_CYCLES);
        test_stat_user_read(PERF_COUNT_HW_BRANCH_MISSES);
        test_stat_user_read(PERF_COUNT_HW_STALLED_CYCLES_FRONTEND);
        test_stat_user_read(PERF_COUNT_HW_STALLED_CYCLES_BACKEND);
//      test_stat_user_read(PERF_COUNT_HW_REF_CPU_CYCLES);
        //test_stat_user_read(PERF_COUNT_HW_REF_CPU_CYCLES);
        //test_stat_read_format();

        __T_END;
        return tests_failed == 0 ? 0 : -1; 
}

}
#endif
//int tests_failed;
//int tests_verbose;
int main(int argc, char** argv)
{
  //#ifdef __cplusplus
  //extern "C"{
  //#endif
  //__T("test evsel", !test_evsel(argc, argv));
  std::ofstream outFile("percentiles.csv");
  outFile << "Set,BG,LS,CB_LEN,50th Percentile,75th Percentile,99th Percentile,99.99th Percentile, Maximum\n";
  parse_args(argc, argv);
  
  benchmarker perf_meas_time("LDPC encoder " + enc_type + " Time", nof_repetitions); //Added by Gayathri
  benchmarker perf_meas_generic("LDPC encoder " + enc_type, nof_repetitions);
  fmt::print("\t-h Show this message TEST  Make\n");
  //uint64_t value_before = read_timer();
  //std::cout << "Value before: " << value_before << std::endl;
  // Capture the current time point
 // std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();
  int iteration_counter = 0; 
  int cycle_count_prints = 0;
  std::vector<int> cd_arr;
  unsigned outst = 0, soutst = 0, inner = 0;
  std::vector<int> cd_arr1;
  std::vector<uint64_t> minValues;
  std::vector<uint64_t> maxValues;
  std::vector<double> avgValues;
  std::vector<double> percentiles50;
  std::vector<double> percentiles75;
  std::vector<double> percentiles99;
  std::vector<double> percentiles999;

  std::vector<int> bg_arr;
  std::vector<int> ls_arr;
  std::vector<int> cb_arr;
  std::vector<int> subset;
  const int subsetSize = 1000;

  for (const ldpc_base_graph_type& bg : {ldpc_base_graph_type::BG1, ldpc_base_graph_type::BG2}) {
	  outst++;
    for (const lifting_size_t& ls : all_lifting_sizes) {
      soutst++;
      std::shared_ptr<ldpc_encoder_factory> encoder_factory = create_ldpc_encoder_factory_sw(enc_type);
      TESTASSERT(encoder_factory);
      std::unique_ptr<ldpc_encoder> encoder = encoder_factory->create();
      TESTASSERT(encoder);

      // Set base-graph message and codeblock lengths.
      unsigned min_cb_length_bg = 24;
      unsigned max_cb_length_bg = 66;//66 was there in srs. i changed it to 68
  //    unsigned std_cb_length_bg = 68;
      unsigned msg_length_bg    = 22;
      if (bg == srsran::ldpc_base_graph_type::BG2) {
        min_cb_length_bg = 12;
        max_cb_length_bg = 50;
        msg_length_bg    = 10;
//	std_cb_length_bg = 52;
      }

      // Compute lifted messages and codeblock lengths.
      unsigned min_cb_length = min_cb_length_bg * ls;
      unsigned max_cb_length = max_cb_length_bg * ls;
      unsigned msg_length    = msg_length_bg * ls;

      for (unsigned cb_length : {min_cb_length, max_cb_length}) {
	inner++;
        // Generate message data.
	iteration_counter++;
        std::vector<uint8_t> data(msg_length);
        std::generate(data.begin(), data.end(), [&]() { return static_cast<uint8_t>(rgen() & 1); });

        // Generate codeblock.
        std::vector<uint8_t> codeblock(cb_length);
	//std::cout<<"cb_length="<<cb_length;
        srsran::codeblock_metadata::tb_common_metadata cfg_enc = {bg, ls};

        fmt::memory_buffer descr_buffer;
        fmt::format_to(descr_buffer, "BG={} LS={:<3} cb_len={}", bg, ls, cb_length);
	// Measure time performance
        perf_meas_time.new_measure(to_string(descr_buffer), data.size(), [&]() {
          encoder->encode(codeblock, data, cfg_enc);
          do_not_optimize(codeblock);
        }); //Gayathri
	// Measure cycle count before the encoder function
        //std::cout << "Cycle count for iteration " << iteration_counter << ": Before encoder" << std::endl;
       // uint64_t start_cycles = readCycleCount();
        //std::cout << "Cycle count for iteration " << iteration_counter << ": Before encoder" << std::endl;
	uint64_t minValue = UINT64_MAX;
	uint64_t maxValue = 0;
	double sum = 0.0;
	//uint64_t a = 0, b = 0;
	perf_meas_generic.new_measure(to_string(descr_buffer), data.size(), [&]() {
	uint64_t start_cycles = readCycleCount(); 
	__T("test evsel",!test_evsel(argc, argv, event_names[eventnum]));
	//__T("test evsel",!test_evsel(argc, argv));
	 uint64_t startCount = get_start_count(global_evsel, &global_counts); 
	 encoder->encode(codeblock, data, cfg_enc);   
	 uint64_t endCount = get_end_count(global_evsel, &global_counts);
	 shutdown_resources(); 
	
	printf("Cycle count: %" PRIu64 "\n", endCount-startCount);
	//std::cout << "a " << ++a << std::endl; 
          #if 0
	//  uint64_t start_cycles = readCycleCount();
          //extern "C"{
         	  
          libperf_init(libperf_print);
	  struct perf_counts_values counts = { .val = 0 };
	  struct perf_thread_map *threads;
    	  struct perf_evsel *evsel;
	  struct perf_event_mmap_page *pc;
	  struct perf_event_attr attr = { 
	  .type  = PERF_TYPE_HARDWARE,
	  .config = PERF_COUNT_HW_INSTRUCTIONS,

	  .config1 = 0x2,     // Request user access 
	  };  
	  int err;

	 // __u64 start;

	  threads = perf_thread_map__new_dummy();
	  if (!threads) {
		fprintf(stderr, "failed to create threads\n");
	//	return -1; 
	  }   

	  perf_thread_map__set_pid(threads, 0, 0); 

	  evsel = perf_evsel__new(&attr);
	  if (!evsel) {
		fprintf(stderr, "failed to create evsel\n");
	//	return -1; 
	  }   

	  err = perf_evsel__open(evsel, NULL, threads);
	  if (err) {
		fprintf(stderr, "failed to open evsel\n");
	//	return -1; 
	  }   

	  err = perf_evsel__mmap(evsel, 0); 
	  if (err) {
		fprintf(stderr, "failed to mmap evsel\n");
	//	return -1; 
	  }   
	  pc = static_cast<perf_event_mmap_page*>(perf_evsel__mmap_base(evsel, 0, 0));

//	  pc = perf_evsel__mmap_base(evsel, 0, 0); 
	  if (!pc) {
		fprintf(stderr, "failed to get mmapped address\n");
	//	return -1; 
	  }   

	  perf_evsel__read(evsel, 0, 0, &counts);
	  uint64_t start = counts.val; 
	// }
          //std::cout << "start value: " << start << std::endl;   
	  encoder->encode(codeblock, data, cfg_enc);
	//extern "C"{  
	perf_evsel__read(evsel, 0, 0, &counts);
	  uint64_t end = counts.val;
	  //std::cout << "end  value: " << end << std::endl;
	printf("Cycle count: %lu\n", end - start);
	perf_evsel__munmap(evsel);
	perf_evsel__close(evsel);
	perf_evsel__delete(evsel);

	perf_thread_map__put(threads);
//}
#endif
         // Measure cycle count after the encoder function
          uint64_t end_cycles = readCycleCount();
          uint64_t cycle_diff = end_cycles - start_cycles;

	  if (cb_length == max_cb_length){
	       //std::cout << "b " << ++b << std::endl;
	  	cd_arr.push_back(cycle_diff);
		cd_arr1.push_back(cycle_diff);
		bg_arr.push_back(static_cast<int>(bg));
		ls_arr.push_back(ls);
		cb_arr.push_back(cb_length);
		minValue = std::min(minValue, cycle_diff);
		maxValue = std::max(maxValue, cycle_diff);
		sum += cycle_diff;
	  }
	  cycle_count_prints++;
	  //std::cout << "Cycle count for iteration " << iteration_counter << ": After encoder (Print " << cycle_count_prints << ")" << std::endl;
          //std::cout << "Cycle count value: " << cycle_diff << std::endl;
	  do_not_optimize(codeblock);
	 return 0; //gayathri
        });
	double average = sum / cd_arr1.size();
	minValues.push_back(minValue);

	maxValues.push_back(maxValue);
	avgValues.push_back(average);
       // for (unsigned i = 0; i < cd_arr.size(); i++){
	//	std::cout << "Cd arr: " << cd_arr[i] << std::endl;}
	//cd_arr.clear();
        cd_arr1.clear();
	// Measure cycle count after the encoder function
       // uint64_t end_cycles = readCycleCount();
       // uint64_t cycle_diff = end_cycles - start_cycles;
        // Print cycle count and iteration information
       // cycle_count_prints++;
       // std::cout << "Cycle count for iteration " << iteration_counter << ": After encoder (Print " << cycle_count_prints << ")" << std::endl;
       // std::cout << "Cycle count value: " << cycle_diff << std::endl;
       //	std::cout << "Cycle count: " << cycle_diff << std::endl;
      }
    }
  }
  std::cout << "size of cd_arr " << cd_arr.size() << std::endl;
  /*for (size_t i = 0; i < minValues.size(); i++) {
    if (i % 2 != 0) {
        std::vector<int> cycleCounts; // Store the cycle counts for this combination
        for (size_t j = i; j < cd_arr.size(); j += 102) {
            cycleCounts.push_back(cd_arr[j]);
        }
        std::cout << "Combination " << i / 2 + 1 << " i=" << i << " j=" << i << " cd_arr.size()=" << cd_arr.size() << std::endl;	
	std::cout << "Combination " << i / 2 + 1 << " Cycle Counts: ";
        for (int cycleCount : cycleCounts) {
            std::cout << cycleCount << " ";
        }
        std::cout << std::endl;

        // Calculate the 50th percentile for this combination
       // double percentile50 = calculatePercentile50(cycleCounts);

        // Store the result in percentiles50
       // percentiles50.push_back(percentile50);
        double percentile50 = calculatePercentile(cycleCounts, 0.5);   // 50th percentile
	double percentile75 = calculatePercentile(cycleCounts, 0.75);  // 75th percentile
	double percentile99 = calculatePercentile(cycleCounts, 0.99);  // 99th percentile
	double percentile999 = calculatePercentile(cycleCounts, 0.999); // 99.9th percentile

        // Print the result
        //std::cout << "Combination " << i / 2 + 1 << ": 50th Percentile = " << percentile50 << std::endl;
	// Print the results
        std::cout << "Combination " << i / 2 + 1 << ": 50th Percentile = " << percentile50 << std::endl;
        std::cout << "Combination " << i / 2 + 1 << ": 75th Percentile = " << percentile75 << std::endl;
        std::cout << "Combination " << i / 2 + 1 << ": 99th Percentile = " << percentile99 << std::endl;
        std::cout << "Combination " << i / 2 + 1 << ": 99.9th Percentile = " << percentile999 << std::endl;
	//cd_arr.clear();
	// Save percentiles to the CSV file
        outFile << "Set " << i << "," << percentile50 << "," << percentile75 << "," << percentile99 << "," << percentile999 << "\n";
    }
  } */
  for (size_t i = 0; i < cd_arr.size(); i++) {
	  subset.push_back(cd_arr[i]);
	  if (subset.size() == subsetSize) {
		//for(unsigned j =0;j<subset.size();j++){
		      //std::cout << subset[j] << " ";}
		//std::cout << std::endl;	
	  	//std::cout << "bg = " << bg_arr[i] << " ls = " << ls_arr[i] << " cb_len = " << cb_arr[i] << std::endl;
		double percentile50 = calculatePercentile(subset, 0.5);   // 50th percentile
         	double percentile75 = calculatePercentile(subset, 0.75);  // 75th percentile
         	double percentile99 = calculatePercentile(subset, 0.99);  // 99th percentile
         	double percentile9999 = calculatePercentile(subset, 0.9999); // 99.99th percentile
         	double percentile100 = calculatePercentile(subset, 1); // 100th percentile
         	std::cout << "Combination " << i / 1000 + 1 << ": 50th Percentile = " << percentile50 << std::endl;
         	std::cout << "Combination " << i / 1000 + 1 << ": 75th Percentile = " << percentile75 << std::endl;
         	std::cout << "Combination " << i / 1000 + 1 << ": 99th Percentile = " << percentile99 << std::endl;
         	std::cout << "Combination " << i / 1000 + 1 << ": 99.99th Percentile = " << percentile9999 << std::endl;
         	std::cout << "Combination " << i / 1000 + 1 << ": 100th Percentile = " << percentile100 << std::endl;
         	// Save percentiles to the CSV file
		outFile << "Set " << i << ","  << bg_arr[i] << "," << ls_arr[i] << "," << cb_arr[i] << "," << percentile50 << "," << percentile75 << "," << percentile99 << "," << percentile9999 <<  "," << percentile100 <<"\n";
		subset.clear();
}
} 

  cd_arr.clear();
  for (unsigned i = 0; i < minValues.size(); i++){
     if(i%2 !=0) 
     	std::cout << "Set " <<i<<": -> Minimum = " << minValues[i]<<", Maximum = "<<maxValues[i]<<", Average = "<<avgValues[i]<< std::endl;
  }

  //std::cout<<"Outermost: "<<outst<<" Second Outermost: " <<soutst<<" Inner: "<<inner<<std::endl;
  //uint64_t value_after = read_timer();
  //std::cout << "Value after: " << value_after << std::endl;
  // Capture another time point after your code execution
//  std::chrono::time_point<std::chrono::steady_clock> end = std::chrono::steady_clock::now();

    // Calculate the elapsed time in nanoseconds
 // auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

    // Print the elapsed time
 // std::cout << "Elapsed time: " << elapsed.count() << " ns" << std::endl;
  // Print time measurements in nanoseconds
  perf_meas_time.print_percentiles_time("nanoseconds"); //Gayathri
  perf_meas_generic.print_percentiles_throughput("bits");
  outFile.close();

}
