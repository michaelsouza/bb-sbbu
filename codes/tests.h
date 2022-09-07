#include "bb.h"
#include "gtest/gtest.h"

TEST( NMR, testA ) {
   NMR nmr( "../DATA_TEST/testA.nmr" );

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

   std::set<int> e1_sid( { 1, 2 } );
   std::set<int> e2_sid( { 2, 3 } );
   std::set<int> e3_sid( { 4 } );

   EXPECT_EQ( E[ 1 ].m_sid, e1_sid );
   EXPECT_EQ( E[ 2 ].m_sid, e2_sid );
   EXPECT_EQ( E[ 3 ].m_sid, e3_sid );
}

TEST( order_sbbu, testA ) {
   NMR nmr( "../DATA_TEST/testA.nmr" );
   std::vector<int> orderSBBU;
   auto costSBBU = order_sbbu( nmr, orderSBBU );
   std::vector<int> order( { 1, 2, 3 } );
   EXPECT_EQ( orderSBBU, order );
}

TEST( order_brute, testA ) {
   NMR nmr( "../DATA_TEST/testA.nmr" );
   std::vector<int> orderOPT;
   auto costOPT = order_brute( nmr, orderOPT );
   EXPECT_EQ( costOPT, 168 );
}

TEST( order_brute, testB ) {
   NMR nmr( "../DATA_TEST/testB.nmr" );
   std::vector<int> orderOPT;
   auto costOPT = order_brute( nmr, orderOPT );
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

void compareBruteWithBB( std::string fnmr ) {
   NMR nmr( fnmr );
   BB bb( nmr );
   std::vector<int> orderOPT;
   auto costOPT = order_brute( nmr, orderOPT );
   auto costBB = bb.solve();
   EXPECT_EQ( costOPT, costBB );
}

TEST( BB, testCostOptimality ) {
   compareBruteWithBB( "../DATA_TEST/testA.nmr" );
   compareBruteWithBB( "../DATA_TEST/testB.nmr" );
   compareBruteWithBB( "../DATA_TEST/testC.nmr" );
   compareBruteWithBB( "../DATA_TEST/testD.nmr" );
   compareBruteWithBB( "../DATA_TEST/testE.nmr" );
}