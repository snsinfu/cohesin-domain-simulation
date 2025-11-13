#include <fstream>
#include <string>
#include <stdexcept>

#include "io.hpp"


std::string
load_text(std::string const& filename)
{
    std::ifstream file(filename);
    std::string text;
    if (!std::getline(file, text, '\0')) {
        throw std::runtime_error("failed to load file " + filename);
    }
    return text;
}
