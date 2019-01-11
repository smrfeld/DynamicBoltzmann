#include <algorithm>
#include <iostream>
#include <cstring>

// for mmap:
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <sstream>
#include <cstdio>
#include <ctime>
#include <map>

const char* map_file(const char* fname, size_t& length);

int main()
{
    clock_t t0 = clock();    

    size_t length;
    auto f = map_file("read_lattice.txt", length);
    auto l = f + length;

    std::string x="",y="",z="";
    std::string sp="";
    std::istringstream iss;

    uintmax_t m_numLines = 0;
    while (f && f!=l)
        if ((f = static_cast<const char*>(memchr(f, '\n', l-f)))) {
            m_numLines++;
            if (!((f != NULL) && (f[0] == '\0'))) {
                iss = std::istringstream(f);
                iss >> x;
                iss >> y;
                iss >> z;
                iss >> sp;
                std::cout << "Line: " << m_numLines << " " << x << " " << y << " " << z << " " << sp << std::endl;

                x = "";
                y = "";
                z = "";
                sp = "";
            };
            f++;
        };

    std::cout << "m_numLines = " << m_numLines << "\n";

    clock_t t1 = clock();    
    std::cout << ( t1 - t0 ) / (double) CLOCKS_PER_SEC << std::endl;

}

void handle_error(const char* msg) {
    perror(msg); 
    exit(255);
}

const char* map_file(const char* fname, size_t& length)
{
    int fd = open(fname, O_RDONLY);
    if (fd == -1)
        handle_error("open");

    // obtain file size
    struct stat sb;
    if (fstat(fd, &sb) == -1)
        handle_error("fstat");

    length = sb.st_size;

    const char* addr = static_cast<const char*>(mmap(NULL, length, PROT_READ, MAP_PRIVATE, fd, 0u));
    if (addr == MAP_FAILED)
        handle_error("mmap");

    // TODO close fd at some point in time, call munmap(...)
    return addr;
}