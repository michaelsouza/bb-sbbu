#include "bb.h"
#include "gtest/gtest.h"

class NMRTest : public ::testing::Test {
protected:
   void SetUp() override {}
};

TEST_F( NMRTest, TestA ) {
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
   for (auto &&kv : E)
   {
      auto eid = kv.first;
      auto& e = kv.second;
      for( auto &&sid : e.m_sid)
      {
         auto &s = S[sid];
         EXPECT_LE(e.m_i + 3, s.m_i);
         EXPECT_LE(e.m_i + 3, s.m_j);
         EXPECT_LE(s.m_i, e.m_j);
         EXPECT_LE(s.m_j, e.m_j);
      }
   }
   
}

TEST( order_sbbu, SolveA ) {
   NMR nmr( "../DATA_TEST/testA.nmr" );
   std::vector<int> orderSBBU;
   auto costSBBU = order_sbbu( nmr, orderSBBU );
   std::vector<int> order( { 1, 2, 3 } );
   EXPECT_EQ( orderSBBU, order );
}

class BBTest : public ::testing::Test {
protected:
   void SetUp() override {}
};

TEST_F( BBTest, SolveA ) {
   NMR nmr( "../DATA_TEST/testA.nmr" );
   BB bb( nmr );
   auto costBB = bb.solve();
   EXPECT_EQ( costBB, 168 );
   int id1, id2;
   for ( int i = 0; i < bb.m_order.size(); i++ ) {
      auto eid = bb.m_order[ i ];
      if ( eid == 1 ) id1 = i;
      if ( eid == 2 ) id2 = i;
   }

   EXPECT_LT( id1, id2 );
}