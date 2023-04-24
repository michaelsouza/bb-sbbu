#include "bb.h"
#include "gtest/gtest.h"

TEST( NMR, testA ) {
   NMR nmr( "../data/nmr_test/testA_chain_A_dmax_5.nmr" );

   std::vector<NMRSegment> segments;
   segments.push_back( NMRSegment( 4, 5 ) );
   segments.push_back( NMRSegment( 6, 10 ) );
   segments.push_back( NMRSegment( 11, 15 ) );
   segments.push_back( NMRSegment( 18, 20 ) );
   EXPECT_EQ( segments.size(), nmr.m_segments.size() );
   for ( size_t i = 0; i < segments.size(); i++ )
      EXPECT_EQ( segments[ i ], nmr.m_segments[ i ] );

   auto& E = nmr.m_E;
   auto& S = nmr.m_S;

   std::vector<int> e1_SID( { 1, 2 } );
   std::vector<int> e2_SID( { 2, 3 } );
   std::vector<int> e3_SID( { 4 } );

   EXPECT_EQ( E[ 1 ].m_SID, e1_SID );
   EXPECT_EQ( E[ 2 ].m_SID, e2_SID );
   EXPECT_EQ( E[ 3 ].m_SID, e3_SID );
}

TEST( sbbuSolve, testA ) {
   NMR nmr( "../data/nmr_test/testA_chain_A_dmax_5.nmr" );
   std::vector<int> orderSBBU;
   auto costSBBU = sbbuSolve( nmr, orderSBBU );
   std::vector<int> order( { 1, 2, 3 } );
   EXPECT_EQ( orderSBBU, order );
}

TEST( sbbuSolve, testF ) {
   NMR nmr( "../data/nmr_test/testF_chain_A_dmax_5.nmr" );
   std::vector<int> orderSBBU;
   auto costSBBU = sbbuSolve( nmr, orderSBBU );
   std::vector<int> order( { 1, 2, 5, 4, 3 } );
   EXPECT_EQ( orderSBBU, order );
}

TEST( bruteSolve, testA ) {
   NMR nmr( "../data/nmr_test/testA_chain_A_dmax_5.nmr" );
   std::vector<int> orderOPT;
   auto costOPT = bruteSolve( nmr, orderOPT );
   EXPECT_EQ( costOPT, 168 );
}

TEST( bruteSolve, testB ) {
   NMR nmr( "../data/nmr_test/testB_chain_A_dmax_5.nmr" );
   std::vector<int> orderOPT;
   auto costOPT = bruteSolve( nmr, orderOPT );
   std::vector<int> order( { 3, 2, 1 } );
   EXPECT_EQ( costOPT, 24 );
   EXPECT_EQ( orderOPT, order );
}

TEST( BST, test ) {
   BST b;
   for ( int i = 1; i <= 10; ++i )
      b.add( i );
   int key = b.minKey( true );
   EXPECT_EQ( key, 1 );
   key = b.minKey( true );
   EXPECT_EQ( key, 2 );
   key = b.minKey( true );
   EXPECT_EQ( key, 3 );
   key = b.minKey( true );
   EXPECT_EQ( key, 4 );
   key = b.minKey( true );
   EXPECT_EQ( key, 5 );
   EXPECT_EQ( b.m_size, 5 );
   b.add( 1 );
   b.add( 2 );
   EXPECT_EQ( b.m_size, 7 );
   key = b.minKeyGT( 2, true );
   EXPECT_EQ( key, 6 );
   EXPECT_EQ( b.m_size, 6 );
   key = b.minKey( true );
   EXPECT_EQ( key, 1 );
   EXPECT_EQ( b.m_size, 5 );
   key = b.minKey( true );
   EXPECT_EQ( key, 2 );
   key = b.minKeyGT( 8, true );
   EXPECT_EQ( key, 9 );
   key = b.minKeyGT( 10, true );
   EXPECT_EQ( key, -1 );
   key = b.minKey( true );
   EXPECT_EQ( key, 7 );
   key = b.minKeyGT( 7, true );
   EXPECT_EQ( key, 8 );
   key = b.minKey( true );
   EXPECT_EQ( key, 10 );
   EXPECT_EQ( b.m_size, 0 );
   key = b.minKey( true );
   EXPECT_EQ( key, -1 );
   EXPECT_EQ( b.m_size, 0 );
}

auto compareBruteWithBB( std::string fnmr ) {
   bool verbose = true; // set true to see the output in debug mode
   NMR nmr( fnmr );
   BB bb( nmr );
   std::vector<int> orderOPT;
   weight_t costOPT = bruteSolve( nmr, orderOPT, 3600, WEIGHT_MAX, verbose );
   if( verbose )
      std::cout << "costOPT = " << costOPT << std::endl;
   
   std::vector<int> orderBB;
   weight_t costBB = 0;
   costBB = bb.solve( costBB, orderBB, 3600, verbose );   
   if( verbose )
      std::cout << "costBB = " << costBB << std::endl;

   printf( "%s\n", fnmr.c_str() );
   EXPECT_EQ( costOPT, costBB );
}

TEST( BB, testCostOptimality ) {
   compareBruteWithBB( "../data/nmr_test/testA_chain_A_dmax_5.nmr" );
   compareBruteWithBB( "../data/nmr_test/testB_chain_A_dmax_5.nmr" );
   compareBruteWithBB( "../data/nmr_test/testC_chain_A_dmax_5.nmr" );
   compareBruteWithBB( "../data/nmr_test/testD_chain_A_dmax_5.nmr" );
   compareBruteWithBB( "../data/nmr_test/testE_chain_A_dmax_5.nmr" );
   compareBruteWithBB( "../data/nmr_test/testF_chain_A_dmax_5.nmr" );
}

auto compareBruteWithPT( std::string fnmr ) {
   NMR nmr( fnmr );
   PT pt( nmr );
   std::vector<int> orderOPT;
   auto costOPT = bruteSolve( nmr, orderOPT );
   auto costPT = pt.solve( 3600, false );
   EXPECT_EQ( costOPT, costPT );
   return pt.m_niters;
}

TEST( PT, testCostOptimality ) {
   compareBruteWithPT( "../data/nmr_test/testA_chain_A_dmax_5.nmr" );
   compareBruteWithPT( "../data/nmr_test/testB_chain_A_dmax_5.nmr" );
   compareBruteWithPT( "../data/nmr_test/testC_chain_A_dmax_5.nmr" );
   compareBruteWithPT( "../data/nmr_test/testD_chain_A_dmax_5.nmr" );
   compareBruteWithPT( "../data/nmr_test/testE_chain_A_dmax_5.nmr" );
   compareBruteWithPT( "../data/nmr_test/testF_chain_A_dmax_5.nmr" );
}

void compareBruteWithGreedy( std::string fnmr ) {
   NMR nmr( fnmr );
   PT pt( nmr );
   std::vector<int> orderOPT;
   auto costBF = bruteSolve( nmr, orderOPT );
   auto costGD = greedySolve( nmr, orderOPT );
   EXPECT_EQ( costBF, costGD );
}

TEST( greedySolve, testCostOptimality ) {
   compareBruteWithGreedy( "../data/nmr_test/testA_chain_A_dmax_5.nmr" );
   compareBruteWithGreedy( "../data/nmr_test/testB_chain_A_dmax_5.nmr" );
   compareBruteWithGreedy( "../data/nmr_test/testC_chain_A_dmax_5.nmr" );
   compareBruteWithGreedy( "../data/nmr_test/testD_chain_A_dmax_5.nmr" );
   compareBruteWithGreedy( "../data/nmr_test/testE_chain_A_dmax_5.nmr" );
   compareBruteWithGreedy( "../data/nmr_test/testF_chain_A_dmax_5.nmr" );
}