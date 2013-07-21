#include <iostream>
#include <sstream>

#define BOOST_TEST_MODULE odeint_mpi
#include <boost/test/unit_test.hpp>

#include <boost/numeric/odeint/external/mpi/mpi.hpp>

using namespace boost::numeric::odeint;

boost::mpi::environment env;

BOOST_AUTO_TEST_SUITE( state_test_suite )

BOOST_AUTO_TEST_CASE( state_test )
{
    boost::mpi::communicator world;

    std::vector<size_t> in_data1, in_data2;
    mpi_state< std::vector<size_t> > state1(world), state2(world);

    // generate data on master
    if(world.rank() == 0) {
        in_data1.resize(31);
        in_data2.resize(33);
        for(size_t i = 0 ; i < in_data2.size() ; i++)
            in_data2[i] = i;
    }

    // copy to nodes
    copy( in_data1, state1 );
    copy( in_data2, state2 );

    {
        std::ostringstream ss;
        ss << "state[" << world.rank() << "] {"
           << state1.data.size() << ", "
           << state2.data.size() << "}\n";
        std::clog << ss.str() << std::flush;
    }

    // compare size
    BOOST_REQUIRE( !same_size( state1, state2 ) );

    // resize state1 to match state2.
    resize( state1, state2 );

    {
        std::ostringstream ss;
        ss << "state[" << world.rank() << "] 1:"
           << state1.data.size() << " 2:"
           << state2.data.size() << "\n";
        std::clog << ss.str() << std::flush;
    }

    // compare size
    BOOST_REQUIRE( same_size( state1, state2 ) );

    // copy state2 to state1
    copy( state2, state1 );

    BOOST_REQUIRE_EQUAL_COLLECTIONS(state1.data.begin(), state1.data.end(),
        state2.data.begin(), state2.data.end());
}


BOOST_AUTO_TEST_SUITE_END()

