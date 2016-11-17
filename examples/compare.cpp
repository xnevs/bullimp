/**
 * Use the VFLib graph database to comapre
 * UllImp to boost's VF2 implementation.
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>

#include <boost/graph/adjacency_list.hpp> 

#include <boost/graph/vf2_sub_graph_iso.hpp>

#include "../graph/ullimp_sub_graph_iso.hpp"


unsigned int read1(std::istream & in) {
    unsigned char a, b;
    a = static_cast<unsigned char>(in.get());
    b = static_cast<unsigned char>(in.get());
    return a | (b << 8);
}

template<typename R>
R read_amalfi(std::istream & in) {
    using namespace boost;

    auto n = read1(in);
    R g{n};
    for(int i=0; i<n; i++) {
        auto cnt = read1(in);
        for(int j=0; j<cnt; j++) {
            auto k = read1(in);
            add_edge(i,k,g);
        }
    }
    return g;
}

using graph_type =  boost::adjacency_list<boost::setS, boost::vecS, boost::bidirectionalS>; 

int main(int argc, char *argv[]) { 
    using namespace std;

    if(argc < 2) {
        cout << "Usage: ./compare test.AXX test.BXX" << endl;
        return 1;
    }

    ifstream in{argv[1], ios::in|ios::binary};
    if(!in.is_open()) return -1;

    auto pattern = read_amalfi<graph_type>(in);
    in.close();

    in.open(argv[2], ios::in|ios::binary);
    if(!in.is_open()) return -1;

    auto graph = read_amalfi<graph_type>(in);
    in.close();

    int count;
    chrono::time_point<chrono::steady_clock> start, end;
    chrono::milliseconds elapsed;

    count = 0;

    start = chrono::steady_clock::now();

    ullimp::ullimp_subgraph_iso(
            pattern,
            graph,
            [&count](auto const &a, auto const &b){
                ++count;
                return true;
            });

    end = chrono::steady_clock::now();
    elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);

    cout << setw(8) << "ullimp: " << count << ", " << "time = " << elapsed.count() << " ms" << endl;


    count = 0;

    start = chrono::steady_clock::now();

    boost::vf2_subgraph_iso(
            pattern,
            graph,
            [&count](auto const &a, auto const &b){
                ++count;
                return true;
            });

    end = chrono::steady_clock::now();
    elapsed = chrono::duration_cast<chrono::milliseconds>(end - start);

    cout << setw(8) << "vf2: " << count << ", " << "time = " << elapsed.count() << " ms" << endl;
}
