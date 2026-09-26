# Boost
find_package(Boost REQUIRED COMPONENTS iostreams)

# OpenMP
find_package(OpenMP REQUIRED)

# Add ISSLCreateIndex executable
add_executable (ISSLCreateIndex ISSLCreateIndex.cpp ISSLCreateIndex.hpp)
target_link_libraries(ISSLCreateIndex PRIVATE Boost::dynamic_linking Boost::iostreams OpenMP::OpenMP_CXX)