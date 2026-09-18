# Boost
find_package(Boost REQUIRED COMPONENTS regex iostreams)

# OpenMP
find_package(OpenMP REQUIRED)

# Add ExtractOfftarget executable
add_executable (ISSLScoreOfftargetsMMF ISSLScoreOfftargetsMMF.cpp ISSLScoreOfftargetsMMF.hpp)
target_compile_options(ISSLScoreOfftargetsMMF PRIVATE -O3 -march=native)
target_link_libraries(ISSLScoreOfftargetsMMF PRIVATE Boost::boost Boost::dynamic_linking Boost::regex Boost::iostreams OpenMP::OpenMP_CXX utils otScorePenalties)