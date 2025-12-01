#include <cstdint>
#include <exception>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>

#include <getopt.hpp>

#include "simulation_config.hpp"
#include "simulation_driver.hpp"

#include "misc/io.hpp"


/** Program mode specified by user. */
enum class program_mode
{
    simulation,
    help,
};


/** Command-line arguments and option flags. */
struct program_options
{
    program_mode                 mode;
    std::string                  config_filename;
    std::optional<std::string>   output_filename;
    std::optional<std::uint64_t> random_seed;
};


static void              show_usage();
static program_options   parse_options(int argc, char** argv);
static simulation_config make_config(program_options const& options);
static simulation_config load_config(std::string const& filename);


int
main(int argc, char** argv)
{
    try {
        program_options const options = parse_options(argc, argv);

        switch (options.mode) {
        case program_mode::simulation:
            simulation_driver(make_config(options)).run();
            break;

        case program_mode::help:
            show_usage();
            break;
        }

        return 0;
    } catch (std::exception const& err) {
        std::cerr << "error: " << err.what() << '\n';
        return 1;
    }
}


void
show_usage()
{
    std::cout <<
        "Chromatin domain simulator\n"
        "usage: main [-osh] <config>\n"
        "\n"
        "  <config>     JSON file specifying simulation parameters\n"
        "\n"
        "options:\n"
        "  -o <output>  override output HDF5 filename ('output_filename' config key)\n"
        "  -s <seed>    override random seed ('random_seed' config key)\n"
        "  -h           print this usage message and exit\n";
}


/** Parses command-line arguments and build a `program_options` structure. */
program_options
parse_options(int argc, char** argv)
{
    program_options options;
    options.mode = program_mode::simulation;

    // Optional arguments
    cxx::getopt getopt;

    for (int opt; (opt = getopt(argc, argv, "o:s:h")) != -1; ) {
        switch (opt) {
        case 'o':
            options.output_filename = getopt.optarg;
            break;

        case 's':
            options.random_seed = std::stoull(getopt.optarg);
            break;

        case 'h':
            options.mode = program_mode::help;
            // Do not parse other arguments if help is requested.
            return options;

        case '?':
            throw std::runtime_error("unrecognized command-line option");
        }
    }

    argc -= getopt.optind;
    argv += getopt.optind;

    // Positional arguments
    if (argc != 1) {
        throw std::runtime_error("missing config file");
    }

    options.config_filename = argv[0];

    return options;
}


/** Builds simulation configuration from given program options. */
simulation_config
make_config(program_options const& options)
{
    simulation_config config = load_config(options.config_filename);

    if (options.output_filename) {
        config.sampling.output_filename = *options.output_filename;
    }

    if (options.random_seed) {
        config.sampling.random_seed = *options.random_seed;
    }

    return config;
}


/** Loads simulation configuration from a JSON file. */
simulation_config
load_config(std::string const& filename)
{
    std::string const text = load_text(filename);
    try {
        return parse_simulation_config(text);
    } catch (std::exception const& err) {
        throw std::runtime_error("failed to parse config file - " + std::string(err.what()));
    }
}
