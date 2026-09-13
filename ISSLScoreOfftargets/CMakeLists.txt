# Boost
find_package(Boost REQUIRED COMPONENTS regex iostreams)

# OpenMP
find_package(OpenMP REQUIRED)

# Add ExtractOfftarget executable
add_executable (ISSLScoreOfftargets ISSLScoreOfftargets.cpp ISSLScoreOfftargets.hpp)
target_compile_options(ISSLScoreOfftargets PRIVATE -O3 -march=native)
target_link_libraries(ISSLScoreOfftargets PRIVATE Boost::boost Boost::dynamic_linking Boost::regex Boost::iostreams OpenMP::OpenMP_CXX utils otScorePenalties)