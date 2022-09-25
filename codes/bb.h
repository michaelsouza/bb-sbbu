#include <algorithm>
#include <chrono>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <regex>
#include <set>
#include <sstream>
#include <stdarg.h>
#include <stdexcept>
#include <string>
#include <vector>

using weight_t = unsigned long long int;
#define WEIGHT_MAX ( ULONG_LONG_MAX / 2 )
#define TIME_NOW() ( std::chrono::steady_clock::now() )
// Elapsed Time in Seconds (ETS)
#define ETS( tic ) ( std::chrono::duration_cast<std::chrono::seconds>( TIME_NOW() - tic ).count() )

class NMRSegment {
private:
   static int SID;

public:
   int m_sid;
   int m_i;
   int m_j;
   weight_t m_weight;
   std::vector<int> m_EID;

   NMRSegment() {
      m_sid = -1;
      m_i = -1;
      m_j = -1;
      m_weight = -1;
   }

   NMRSegment( const NMRSegment& s )
       : m_EID( s.m_EID ) {
      m_sid = s.m_sid;
      m_i = s.m_i;
      m_j = s.m_j;
      m_weight = s.m_weight;
   }

   NMRSegment( int i, int j ) {
      m_sid = ++NMRSegment::SID;
      m_i = i;
      m_j = j;
      updateWeight();
   }

   static void resetSID() {
      NMRSegment::SID = 0;
   }

   void addEid( int eid ) {
      // ordered insertion
      auto it = std::upper_bound( m_EID.begin(), m_EID.end(), eid );
      // element not found
      if ( it == m_EID.end() )
         m_EID.insert( it, eid );
   }

   void updateWeight() {
      weight_t p = m_j - m_i + 1;
      if ( p > 63 ) throw std::invalid_argument( "Segment weight is too large to be represented (overflow)." );
      m_weight = 1 << p;
   }

   bool operator==( NMRSegment const& rhs ) const {
      return m_i == rhs.m_i && m_j == rhs.m_j;
   }
};

int NMRSegment::SID;

class NMREdge {
private:
   static int EID; // class static variable

public:
   int m_eid;
   int m_i;
   int m_j;
   std::vector<int> m_SID;

   NMREdge() {
      m_eid = -1;
      m_i = -1;
      m_j = -1;
   }

   NMREdge( int i, int j ) {
      m_eid = ++NMREdge::EID; // edge id
      m_i = i;
      m_j = j;
   }

   void addSid( int sid ) {
      // ordered insertion
      auto it = std::upper_bound( m_SID.begin(), m_SID.end(), sid );
      // element not found
      if ( it == m_SID.end() )
         m_SID.insert( it, sid );
   }

   bool checkCover( NMRSegment& s ) {
      return ( m_i + 3 <= s.m_i ) && ( s.m_j <= m_j );
   }

   static void resetEID() {
      NMREdge::EID = 0;
   }
};

int NMREdge::EID;

class NMR {
public:
   std::string m_fnmr;
   unsigned int m_nnodes;
   std::vector<NMREdge> m_edges;
   std::vector<NMREdge> m_pruneEdges;
   std::vector<NMRSegment> m_segments;
   std::map<int, NMREdge> m_E;
   std::map<int, NMRSegment> m_S;

   NMR( std::string fnmr ) {
      m_fnmr = fnmr;

      NMREdge::resetEID(); // set EID=0

      // open fnmr file
      std::ifstream fid( fnmr.c_str() );
      if ( !fid.good() )
         throw std::runtime_error( "Could not open fnmr file." );

      // read fnmr file row by row
      int i, j;
      std::string row;
      while ( std::getline( fid, row ) ) {
         std::stringstream ss( row );
         ss >> i;
         ss >> j;
         m_edges.push_back( NMREdge( i, j ) );
      }

      if ( m_edges.size() == 0 )
         throw std::runtime_error( "No edges found." );

      m_nnodes = 0;
      for ( auto&& e : m_edges ) {
         if ( m_nnodes < e.m_j ) m_nnodes = e.m_j;
         if ( e.m_j > e.m_i + 3 ) m_pruneEdges.push_back( e );
      }
      setSegments();
      setOrderingData();
   }

private:
   void setSegments() {
      NMRSegment::resetSID();
      // I: sorted list of all atoms covered by prune edges
      std::vector<int> I;
      for ( auto&& e : m_pruneEdges )
         for ( int i = e.m_i + 3; i <= e.m_j; ++i )
            I.push_back( i );
      std::sort( I.begin(), I.end() );
      I.erase( std::unique( I.begin(), I.end() ), I.end() );

      // E[i]: set of the eid's of the edges covering atom 'i'
      std::map<int, std::set<int>> E;
      for ( auto&& i : I ) E[ i ] = {};

      for ( auto&& e : m_pruneEdges )
         for ( int i = e.m_i + 3; i <= e.m_j; i++ )
            E[ i ].insert( e.m_eid );

      // Consecutive atoms covered by the same set of edges belong to the same segment
      {
         int s_i = I[ 0 ], s_j = I[ 0 ];
         for ( auto&& j : I ) {
            if ( E[ s_i ] == E[ j ] )
               s_j = j;
            else {
               m_segments.push_back( NMRSegment( s_i, s_j ) );
               s_i = s_j = j;
            }
         }
         m_segments.push_back( NMRSegment( s_i, s_j ) );
      }

      // O(len(S) * len(m_pruneEdges))
      for ( auto&& s : m_segments )
         for ( auto&& e : m_pruneEdges )
            if ( e.checkCover( s ) ) {
               e.addSid( s.m_sid );
               s.addEid( e.m_eid );
            }
   }

   void setOrderingData() {
      for ( auto&& e : m_pruneEdges )
         m_E[ e.m_eid ] = e;
      for ( auto&& s : m_segments )
         m_S[ s.m_sid ] = s;
   }
};

weight_t costOrder( std::vector<int>& order, std::map<int, NMREdge>& E, std::map<int, NMRSegment>& S, weight_t costUB = WEIGHT_MAX ) {
   weight_t total_cost = 0; // total cost
   // sid is added to B by the first edge that covers it.
   std::set<int> B;
   for ( auto&& eid : order ) {
      weight_t edge_cost = 1;
      const auto& e = E[ eid ];
      for ( auto&& sid : e.m_SID ) {
         // sid not in S, it does not need to be covered
         if ( S.find( sid ) == S.end() ) continue;
         // sid already covered
         if ( B.find( sid ) != B.end() ) continue;
         // covering sid
         const auto& s = S[ sid ];
         edge_cost *= s.m_weight;
         B.insert( sid );
      }
      // the cost of an edge that covers no segment is zero (not one)
      total_cost += edge_cost > 1 ? edge_cost : 0;
      if ( total_cost >= costUB )
         return WEIGHT_MAX;
   }
   return total_cost;
}

weight_t sbbuSolve( NMR& nmr, std::vector<int>& order ) {
   auto& E = nmr.m_E;

   order.clear();
   for ( auto&& kv : E )
      order.push_back( kv.first ); // list of edges eid

   std::sort( order.begin(), order.end(),
       [ &E ]( const int& eidA, const int& eidB ) -> bool {
          const auto& eA = E[ eidA ];
          const auto& eB = E[ eidB ];
          if ( eA.m_j == eB.m_j )
             return eA.m_i > eB.m_i;
          return eA.m_j < eB.m_j;
       } );

   return costOrder( order, E, nmr.m_S );
}

weight_t bruteSolve( NMR& nmr, std::vector<int>& orderOPT ) {
   auto& E = nmr.m_E;
   auto& S = nmr.m_S;
   weight_t costOPT = WEIGHT_MAX;

   std::vector<int> order;
   for ( auto&& kv : E )
      order.push_back( kv.first );

   orderOPT.resize( order.size() );

   do {
      auto c = costOrder( order, E, S );
      if ( c < costOPT ) {
         costOPT = c;
         std::copy( order.begin(), order.end(), orderOPT.begin() );
      }
   } while ( std::next_permutation( order.begin(), order.end() ) );

   return costOPT;
}

weight_t costRelax( std::set<int>& U, std::map<int, NMRSegment>& S ) {
   weight_t costTotal = 0;
   for ( auto&& sid : U ) costTotal += S[ sid ].m_weight;
   return costTotal;
}

weight_t costRelax( std::map<int, NMRSegment>& S ) {
   weight_t costTotal = 0;
   for ( auto&& kv : S )
      costTotal += S[ kv.first ].m_weight;
   return costTotal;
}

/**
 * @brief Binary Search Tree (BST)
 *
 */
class BSTNode {
public:
   int m_key;
   BSTNode* m_lft; // left child
   BSTNode* m_rht; // right child
   BSTNode* m_p;   // parent

   BSTNode( int key = -1 ) {
      m_key = key;
      m_lft = nullptr;
      m_rht = nullptr;
      m_p = nullptr;
   }
};

class BST {
public:
   BSTNode* m_root;
   size_t m_size;
   BST() {
      m_root = nullptr;
      m_size = 0;
   }

   void inorderTreeWalk( BSTNode* x = nullptr ) {
      if ( x == nullptr ) { // base case
         inorderTreeWalk( m_root );
         return;
      }
      if ( x->m_lft ) inorderTreeWalk( x->m_lft );
      printf( "%d, ", x->m_key );
      if ( x->m_rht ) inorderTreeWalk( x->m_rht );
   }

   BSTNode* treeSearch( int key ) {
      BSTNode* x = m_root;
      while ( x != nullptr && key != x->m_key )
         x = key < x->m_key ? x->m_lft : x->m_rht;
      return x;
   }

   BSTNode* minNode( BSTNode* x = nullptr ) {
      if ( m_size == 0 ) return nullptr;
      if ( x == nullptr ) x = m_root;
      while ( x->m_lft != nullptr ) x = x->m_lft;
      return x;
   }

   int minKey( bool deleteNode = false ) {
      auto x = minNode();
      if ( x == nullptr ) return -1;
      auto key = x->m_key;
      if ( deleteNode ) rem( x, true );
      return key;
   }

   BSTNode* maxNode( BSTNode* x = nullptr ) {
      if ( m_size == 0 ) return nullptr;
      if ( x == nullptr ) x = m_root;
      while ( x->m_rht != nullptr ) x = x->m_rht;
      return x;
   }

   int maxKey( bool deleteNode = false ) {
      auto x = maxNode();
      if ( x == nullptr ) return -1;
      auto key = x->m_key;
      if ( deleteNode ) rem( x, true );
      return key;
   }

   BSTNode* minNodeGT( int key ) {
      BSTNode* z = nullptr;
      BSTNode* x = m_root;
      int minKey = INT_MAX;
      while ( x != nullptr ) {
         if ( x->m_key > key ) {
            minKey = key;
            z = x;
            x = x->m_lft;
         }
         else { // x->m_key >= key
            x = x->m_rht;
         }
      }
      return z;
   }

   int minKeyGT( int key, bool deleteNode = false ) {
      auto x = minNodeGT( key );
      // return an invalid key
      if ( x == nullptr ) return -1;
      key = x->m_key;
      if ( deleteNode ) rem( x, true );
      return key;
   }

   void add( int key ) {
      ++m_size;
      BSTNode* x = m_root;
      BSTNode* y = nullptr; // will be the parent of key
      while ( x != nullptr ) {
         y = x;
         x = x->m_key > key ? x->m_lft : x->m_rht;
      }
      BSTNode* z = new BSTNode( key );
      z->m_p = y;
      if ( y == nullptr )
         m_root = z; // tree was empty
      else if ( key < y->m_key )
         y->m_lft = z;
      else
         y->m_rht = z;
   }

   void transplant( BSTNode* u, BSTNode* v ) {
      // update u's parent
      if ( u->m_p == nullptr ) {
         m_root = v;
      }
      else if ( u == u->m_p->m_lft ) {
         u->m_p->m_lft = v;
      }
      else {
         u->m_p->m_rht = v;
      }

      // update v's parent
      if ( v != nullptr )
         v->m_p = u->m_p;
   }

   void rem( BSTNode* z, bool deleteNode = false ) {
      if ( z == nullptr ) return;
      --m_size;
      if ( z->m_lft == nullptr ) { // (possibly) only rht child
         transplant( z, z->m_rht );
      }
      else if ( z->m_rht == nullptr ) { // only lft child
         transplant( z, z->m_lft );
      }
      else {
         auto y = minNode( z->m_rht );
         if ( y != z->m_rht ) {
            transplant( y, y->m_rht );
            y->m_rht = z->m_rht;
            y->m_rht->m_p = y;
         }
         transplant( z, y );
         y->m_lft = z->m_lft;
         y->m_lft->m_p = y;
      }
      if ( deleteNode )
         delete z;
   }
};

class BBP {
public:
   // current order
   std::vector<int> m_ord;
   // index of the element last element inserted
   int m_idx;
   // n:normal, p:prune
   char m_state;
   // binary search tree of available eid's
   BST m_bst;
   // m_nUncSID[eid]: number of uncovered sid's associated to eid
   std::vector<int> m_nUncSID;
   // m_nCovEID[sid]: number of eid's in the order that cover sid
   std::vector<int> m_nCovEID;
   // m_binBst[eid]: true iff eid was inserted in m_t
   std::vector<bool> m_binBst;
   // m_binOrd[eid]: true iff eid was inserted in m_o
   std::vector<bool> m_binOrd;
   std::map<int, NMREdge>& m_E;
   std::map<int, NMRSegment>& m_S;

   BBP( std::map<int, NMREdge>& E, std::map<int, NMRSegment>& S )
       : m_E( E ), m_S( S ) {
      m_state = 'n';
      m_idx = -1;

      m_ord.resize( E.size() );
      std::fill( m_ord.begin(), m_ord.end(), -1 );

      int eidMax = 0;
      for ( auto&& kv : E ) {
         m_bst.add( kv.first );
         if ( eidMax < kv.first ) eidMax = kv.first;
      }

      int sidMax = 0;
      for ( auto&& kv : S )
         if ( sidMax < kv.first ) sidMax = kv.first;

      m_nUncSID.resize( eidMax + 1 );
      for ( auto&& kv : E )
         m_nUncSID[ kv.first ] = kv.second.m_SID.size();

      m_nCovEID.resize( sidMax + 1 );
      std::fill( m_nCovEID.begin(), m_nCovEID.end(), 0 );

      m_binBst.resize( m_nUncSID.size() ); // key can be one-based
      std::fill( m_binBst.begin(), m_binBst.end(), true );

      m_binOrd.resize( m_binBst.size() );
      std::fill( m_binOrd.begin(), m_binOrd.end(), false );
   }

   inline int eidMin() {
      auto key = m_bst.minKey( true );
      m_binBst[ key ] = false;
      return key;
   }

   inline int eidMinGT( int eid ) {
      eid = m_bst.minKeyGT( eid, true );
      if ( eid > 0 ) m_binBst[ eid ] = false;
      return eid;
   }

   inline void addBst( int eid ) {
      if ( eid > 0 && m_binOrd[ eid ] )
         throw std::runtime_error( "Trying to add to m_bst an eid already in m_ord." );
      if ( m_binBst[ eid ] )
         throw std::runtime_error( "Trying to add to m_bst an eid already included." );
      m_bst.add( eid );
      m_binBst[ eid ] = true;
   }

   inline void addOrd( const int eid ) {
      if ( eid < 1 ) return;
      if ( m_binBst[ eid ] )
         throw std::runtime_error( "Trying to add to m_ord an eid already in m_bst." );
      if ( m_binOrd[ eid ] )
         throw std::runtime_error( "Trying to add to m_ord an eid already included." );

      m_ord[ ++m_idx ] = eid;
      m_binOrd[ eid ] = true;

      const auto& e = m_E[ eid ];
      for ( auto&& sid : e.m_SID ) {
         const auto& s = m_S[ sid ];
         ++m_nCovEID[ sid ];
         if ( m_nCovEID[ sid ] == 1 )
            for ( auto&& eid : s.m_EID ) --m_nUncSID[ eid ];
      }
   }

   inline int remOrd() {
      const auto eid = m_ord[ m_idx ];
      // set invalid value
      m_ord[ m_idx-- ] = -1;
      if ( eid < 1 ) return eid;

      m_binOrd[ eid ] = false;
      const auto e = m_E[ eid ];
      for ( auto&& sid : e.m_SID ) {
         const auto s = m_S[ sid ];
         --m_nCovEID[ sid ];
         // there is no eid in 'ord' covering sid
         // add all eid associated to it in the 'bst'
         if ( m_nCovEID[ sid ] == 0 ) {
            for ( auto eid : s.m_EID ) {
               ++m_nUncSID[ eid ];
               if ( m_binBst[ eid ] == false && m_binOrd[ eid ] == false )
                  addBst( eid );
            }
         }
      }
      return eid;
   }

   /**
    * @brief Return the index of the last component in the current order. If there is no
    *        next component, the index -1 is returned.
    *
    * @return int
    */
   int next( bool verbose = false ) {
      if ( m_state == 'n' ) {
         const auto eid = eidMin();
         // m_t is empty
         if ( eid == -1 ) {
            m_state = 'p';
            return next();
         }

         // just drop key from m_t and call next
         if ( m_nUncSID[ eid ] == 0 ) return next();

         addOrd( eid );
         return eid;
      }

      if ( m_idx == -1 ) return -1;

      const auto eidOld = remOrd();
      const auto eid = eidMinGT( eidOld );
      if ( eid < 1 || m_nUncSID[ eid ] == 0 ) {
         return next();
      }
      addOrd( eid );
      m_state = 'n';
      return eid;
   }

   void prune() {
      m_state = 'p';
   }
};

class BB {
public:
   unsigned long long int m_niters;
   int m_idx;
   NMR& m_nmr;
   BBP m_p;
   bool m_timeout;
   size_t m_nedges;
   std::vector<int> m_ord;
   std::map<int, NMREdge>& m_E;
   std::map<int, NMRSegment>& m_S;
   // m_c[i]: cost added by the i-th eid of the current order
   std::vector<weight_t> m_c;
   // m_q[eid]: number of uncovered sid's associated to eid

   BB( NMR& nmr )
       : m_nmr( nmr ), m_E( nmr.m_E ), m_S( nmr.m_S ), m_p( nmr.m_E, nmr.m_S ) {
      m_idx = -1;
      m_nedges = m_E.size();
      m_ord.resize( m_nedges );
      m_timeout = false;
      // init c *****************
      m_c.resize( m_nedges );
      std::fill( m_c.begin(), m_c.end(), 0 );
   }

   /**
    * @brief
    *
    * @param C C[id] is the number of edges on m_ord covering the segment sid.
    * @param U set of available segments
    * @return weight_t
    */
   weight_t remOrd( std::vector<int>& C, weight_t& costRLX ) {
      // total cost of removed eid's
      weight_t costTotal = 0;
      while ( m_idx >= m_p.m_idx && m_p.m_ord[ m_idx ] != m_ord[ m_idx ] ) {
         auto eid = m_ord[ m_idx ];
         // set invalid value
         m_ord[ m_idx ] = -1;
         costTotal += m_c[ m_idx ];
         const auto& e = m_E[ eid ];
         for ( auto&& sid : e.m_SID ) {
            --C[ sid ];
            if ( C[ sid ] == 0 ) {
               const auto s = m_S[ sid ];
               // increment the relaxed cost
               costRLX += s.m_weight;
            }
         }
         --m_idx;
      }
      return costTotal;
   }

   /**
    * @brief Add edge eid to the m_order and update C and U.
    *
    * @param eid
    * @param C C[sid] : number of edges on m_order covering segment sid
    * @param U set of the uncovered segments.
    * @return weight_t
    */
   weight_t addOrd( int eid, std::vector<int>& C, weight_t& costRLX,
       const weight_t costACC, const weight_t costUB ) {
      m_ord[ ++m_idx ] = eid;
      weight_t costEid = 1;
      const auto& e = m_E[ eid ];

      for ( auto&& sid : e.m_SID ) {
         ++C[ sid ];
         // the current eid is the only one covering the sid
         if ( C[ sid ] == 1 ) {
            const auto& s = m_S[ sid ];
            // reduce the relaxed cost
            costRLX -= s.m_weight;
            // avoid overflow by not incrementing costEid,
            // if it's too large. the costEid will not be correctly
            // calculated, but the eid will be pruned any way.
            if ( costUB >= costACC + costEid )
               costEid *= s.m_weight;
         }
      }
      costEid = costEid > 1 ? costEid : 0;
      m_c[ m_idx ] = costEid;
      return costEid;
   }

   weight_t solve( size_t tmax = 3600, bool verbose = false ) {
      if ( verbose ) printf( "\n\nsolving %s\n", m_nmr.m_fnmr.c_str() );
      m_niters = 0;
      auto tic = TIME_NOW();
      weight_t costUB = WEIGHT_MAX;
      std::vector<int> orderOPT;
      costUB = sbbuSolve( m_nmr, orderOPT );

      // init m_order
      std::fill( m_ord.begin(), m_ord.end(), -1 );

      // C[sid] : number of edges already included in the order that cover segment sid
      std::vector<int> C( m_S.size() + 1 ); // the segments are one-based
      std::fill( C.begin(), C.end(), 0 );

      // first cost_relax
      auto costRELAX = costRelax( m_S );

      // solution found
      if ( costRELAX == costUB ) {
         std::copy( orderOPT.begin(), orderOPT.end(), m_ord.begin() );
         return costUB;
      }

      // costACC: cost accumulated to the current idx
      weight_t costACC = 0;
      // costRLX: LB of the cost to cover the remaning segments
      // When costRLX == 0, there is no remaining segments to be covered
      weight_t costRLX = costRELAX;
      // costLB = costACC + costRLX
      weight_t costLB = costACC + costRLX;
      // eid_cost: cost of the last edge added
      weight_t costEID = costACC + costRLX;

      int eid = m_p.next();
      // loop through all permutations
      while ( eid > -1 ) {
         ++m_niters;
         auto toc = ETS( tic );
         if ( toc > tmax ) {
            m_timeout = true;
            printf( "> timeoutBB %ld seconds\n", toc );
            break;
         }

         costACC -= remOrd( C, costRLX );
         costEID = addOrd( eid, C, costRLX, costACC, costUB );

         costACC += costEID;
         costLB = costACC + costRLX;

         if ( verbose ) {
            printf( "UB:%8llu, LB:%8llu, CE:%8llu, AC:%8llu, RL:%8llu, ", costUB, costLB, costEID, costACC, costRLX );
            printf( " o:[" );
            for ( int i = 0; i <= m_idx; ++i )
               printf( "%d, ", m_ord[ i ] );
            printf( "]\n" );
         }

         // sanity check (something went wrong)
         if ( costLB < costRELAX )
            throw std::runtime_error( "Impossible lower bound (costLB < costRELAX)." );

         // update solution
         if ( ( costRLX == 0 ) && ( costLB < costUB ) ) {
            costUB = costLB;
            std::copy( m_ord.begin(), m_ord.end(), orderOPT.begin() );
            // best possible solution found
            if ( costRELAX == costLB ) break;
         }

         // prune when
         // 1) costLB > costUB: the prefix will not generate a improved solution
         // 2) costRLX == 0: there is no remaining segment to be solved
         // 3) costEID == 0: push the solved eid to the end of the permutation
         if ( ( costUB <= costLB ) || ( costRLX == 0 ) || ( costEID == 0 ) )
            m_p.prune();

         eid = m_p.next();
      }

      std::copy( orderOPT.begin(), orderOPT.end(), m_ord.begin() );
      return costUB;
   }
};

bool exists( std::string fname ) {
   std::ifstream fid( fname.c_str() );
   return fid.good();
}

void write_log( FILE* fid, const char* fmt, ... ) {
   va_list args;
   va_start( args, fmt );
   vfprintf( fid, fmt, args );
   va_end( args );

   va_start( args, fmt );
   vprintf( fmt, args );
   va_end( args );

   fflush( fid );
}

class PT {
public:
   std::vector<std::set<int>> m_preds;
   std::vector<int> m_ordS;
   NMR& m_nmr;
   int m_niters;
   std::map<int, NMRSegment>& m_S;
   bool m_timeout;
   std::vector<bool> m_b;
   std::vector<int> m_p;
   std::vector<int> m_Ek;
   PT( NMR& nmr )
       : m_nmr( nmr ), m_S( nmr.m_S ), m_b( nmr.m_nnodes, false ) {
      
      init_ordS( );
      
      // Ek[i]: num of uncovered segments of edge 'i'
      m_Ek.resize( nmr.m_edges.size() + 1 );
      for ( auto&& kv : nmr.m_E ) {
         m_Ek[ kv.second.m_eid ] = kv.second.m_SID.size();
      }

      m_preds.resize( m_nmr.m_edges.size() + 1 );
   }

   void init_ordS( ) {
      std::vector<NMREdge> E;
      for ( auto& kv : m_nmr.m_E ) E.push_back( kv.second );
      // ascend sort of edges by degree
      std::sort( E.begin(), E.end(), []( const NMREdge& a, const NMREdge& b ) -> bool {
         return a.m_SID.size() > b.m_SID.size();
      } );
      for ( auto& kv : m_nmr.m_E ) {
         auto& e = kv.second;
         std::vector<NMRSegment*> S;
         // ascend sort of edge segments by degree
         for ( auto&& sid : e.m_SID ) S.push_back( &m_nmr.m_S[ sid ] );
         sort( S.begin(), S.end(), []( const NMRSegment* a, const NMRSegment* b ) -> bool {
            return a->m_EID.size() > b->m_EID.size();
         } );
         for ( auto&& s : S )
            if ( m_b[ s->m_sid ] == false ) {
               m_b[ s->m_sid ] = true;
               m_ordS.push_back( s->m_sid );
            }
      }
      // restore m_b
      for ( auto& sid : m_ordS ) m_b[ sid ] = false;
   }

   void predecessors( int eidA, std::vector<int>& p ) {
      p.clear();
      int i = 0;
      for ( auto& eidB : m_preds[ eidA ] ) {
         p.push_back( eidB );
         m_b[ eidB ] = true;
      }
      while ( i < p.size() ) {
         for ( auto& eidB : m_preds[ p[ i ] ] )
            if ( m_b[ eidB ] == false ) {
               m_b[ eidB ] = true;
               p.push_back( eidB );
            }
         ++i;
      }
      sort( p.begin(), p.end() );
      // restore m_b
      for ( auto& eidB : p ) m_b[ eidB ] = false;
   }

   void available_edges( int sid, std::vector<int>& E ) {
      E.clear();
      auto& EID = m_S[ sid ].m_EID;
      
      if ( EID.size() == 1 ) {
         E.push_back( EID[ 0 ] );
         return;
      }

      for ( size_t i = 0; i < EID.size(); i++ ) {
         const auto& eidA = EID[ i ];
         if ( m_preds[ eidA ].size() == 0 ) {
            E.push_back( eidA );
            continue;
         }
         predecessors( eidA, m_p );
         bool avail = true;
         for ( size_t j = 0; j < EID.size(); j++ ) {
            const auto& eidB = EID[ j ];
            if ( eidA == eidB ) continue;
            const auto it = std::lower_bound( m_p.begin(), m_p.end(), eidB );
            if ( ( it < m_p.end() ) && ( *it ) == eidB ) {
               avail = false;
               break;
            }
         }
         if ( avail ) E.push_back( eidA );
      }
   }

   void add_precedence( int eidA, std::vector<int>& E, std::vector<std::pair<int, int>>& P ) {
      for ( auto&& eidB : E ) {
         if ( eidA == eidB ) continue;
         P.push_back( std::make_pair( eidB, eidA ) );
         m_preds[ eidB ].insert( eidA );
      }
   }

   weight_t edge_cost( std::vector<int>& c_eid, int eid, weight_t costUB ) {
      weight_t cost = 1;
      for ( auto&& sid : m_nmr.m_E[ eid ].m_SID ) {
         if ( c_eid[ sid ] == eid ) {
            const auto& s = m_S[ sid ];
            cost *= s.m_weight;
            if ( cost >= costUB ) {
               cost = costUB;
               break;
            }
         }
      }
      return cost == 1 ? 0 : cost;
   }

   weight_t add_cost( int sid, std::vector<int>& c_eid, weight_t costUB ) {
      weight_t cost = 0;
      for ( auto& eid : m_S[ sid ].m_EID ) {
         if ( m_Ek[ eid ] == 0 ) continue;
         m_Ek[ eid ] -= 1;
         if ( m_Ek[ eid ] == 0 && cost < costUB ) {
            cost += edge_cost( c_eid, eid, costUB );
         }
      }
      return cost;
   }

   void rem_precedence( std::vector<std::pair<int, int>>& P ) {
      for ( const auto& p : P )
         m_preds[ p.first ].erase( p.second );
      P.clear();
   }

   weight_t rem_cost( int level, int sid, std::vector<weight_t>& costADD ) {
      auto costREM = costADD[ level ];
      costADD[ level ] = 0;
      for ( const auto& eid : m_S[ sid ].m_EID )
         ++m_Ek[ eid ];
      return costREM;
   }

   void backtracking( int& level, std::vector<std::vector<int>>& E, std::vector<std::vector<std::pair<int, int>>>& P, std::vector<int>& c_idx, std::vector<int>& c_eid, weight_t& cost, std::vector<weight_t>& costADD ) {
      while ( level >= 0 ) {
         auto sid = m_ordS[ level ];
         cost -= rem_cost( level, sid, costADD );
         c_eid[ sid ] = -1;
         if ( P[ level ].size() > 0 ) rem_precedence( P[ level ] );
         if ( c_idx[ level ] < E[ level ].size() - 1 ) {
            ++c_idx[ level ];
            return;
         }
         E[ level ].clear();
         c_idx[ level ] = 0;
         --level;
      }
   }

   weight_t solve( size_t tmax = 3600, bool verbose = false ) {
      if ( verbose ) printf( "\n\nsolving %s\n", m_nmr.m_fnmr.c_str() );
      m_niters = 0;
      m_timeout = false;
      auto tic = TIME_NOW();
      std::vector<int> orderOPT;
      weight_t costOPT = sbbuSolve( m_nmr, orderOPT );
      auto costRELAX = costRelax( m_S );
      // solution found
      if ( costRELAX == costOPT ) return costOPT;
      // c_eid[sid]: eid covering sid
      std::vector<int> c_eid( m_nmr.m_segments.size() + 1 );
      std::vector<int> c_idx( m_ordS.size(), 0 );
      std::vector<weight_t> costADD( m_ordS.size(), 0 );
      int level = 0, sid, eid;
      weight_t cost = 0;
      std::vector<std::vector<int>> E( m_ordS.size() );
      std::vector<std::vector<std::pair<int, int>>> P( m_ordS.size() );
      auto toc = ETS( tic );
      while ( level >= 0 ) {
         ++m_niters;
         toc = ETS( tic );
         if ( toc > tmax ) {
            m_timeout = true;
            printf( "> timeoutBB %ld seconds\n", toc );
            break;
         }
         sid = m_ordS[ level ];
         if ( E[ level ].size() == 0 ) available_edges( sid, E[ level ] );
         eid = E[ level ][ c_idx[ level ] ];
         c_eid[ sid ] = eid;
         if ( E[ level ].size() >= 2 ) add_precedence( eid, E[ level ], P[ level ] );
         costADD[ level ] = add_cost( sid, c_eid, costOPT );
         cost += costADD[ level ];
         // solution found
         if ( cost < costOPT && level == ( m_ordS.size() - 1 ) ) costOPT = cost;
         // next
         if ( cost < costOPT && level < ( m_ordS.size() - 1 ) )
            ++level;
         else
            backtracking( level, E, P, c_idx, c_eid, cost, costADD );
      }
      return costOPT;
   }
};

int call_solvers( int argc, char* argv[] ) {
   std::string fnmr = "/home/michael/gitrepos/bb-sbbu/DATA_TEST/testC.nmr";
   size_t tmax = 3600;
   bool clean_log = false;
   bool verbose = false;
   for ( int i = 0; i < argc; ++i ) {
      auto arg = argv[ i ];
      if ( strcmp( arg, "-fnmr" ) == 0 ) fnmr = argv[ i + 1 ];
      if ( strcmp( arg, "-tmax" ) == 0 ) tmax = atoi( argv[ i + 1 ] );
      if ( strcmp( arg, "-clean_log" ) == 0 ) clean_log = true;
      if ( strcmp( arg, "-verbose" ) == 0 ) verbose = true;
   }

   auto flog = std::regex_replace( fnmr, std::regex( ".nmr" ), ".log" );
   // check if already has a log file
   if ( !clean_log && exists( flog ) ) {
      printf( "> skip (already solved) %s\n", fnmr.c_str() );
      return EXIT_SUCCESS;
   }
   // create log file
   FILE* fid = fopen( flog.c_str(), "w" );
   if ( fid == nullptr )
      throw std::runtime_error( "Could create the log file." );

   write_log( fid, "> fnmr %s\n", fnmr.c_str() );

   // read instance
   NMR nmr( fnmr );
   auto& E = nmr.m_E;
   auto& S = nmr.m_S;

   write_log( fid, "> verbose ........... %d\n", verbose ? 1 : 0 );
   write_log( fid, "> clean_log ......... %d\n", clean_log ? 1 : 0 );
   write_log( fid, "> tmax (secs) ....... %ld\n", tmax );
   write_log( fid, "> nnodes ............ %d\n", nmr.m_nnodes );
   write_log( fid, "> lenE .............. %d\n", E.size() );
   write_log( fid, "> lenS .............. %d\n", S.size() );

   std::set<int> U;
   for ( auto&& kv : S ) U.insert( kv.first );

   auto costRELAX = costRelax( U, S );
   write_log( fid, "> costRELAX ......... %d\n", costRELAX );

   // call order_sbbu
   std::vector<int> orderSBBU;
   auto tic = TIME_NOW();
   auto costSBBU = sbbuSolve( nmr, orderSBBU );
   auto toc = ETS( tic );
   write_log( fid, "> costSBBU .......... %d\n", costSBBU );
   write_log( fid, "> timeSBBU (secs) ... %ld\n", toc );

   // call order_bb
   BB bb( nmr );
   tic = TIME_NOW();
   auto costBB = bb.solve( tmax, verbose );
   toc = ETS( tic );
   write_log( fid, "> timeoutBB ......... %d\n", bb.m_timeout ? 1 : 0 );
   write_log( fid, "> costBB ............ %d\n", costBB );
   write_log( fid, "> timeBB (secs) ..... %ld\n", toc );

   // call order_bb if needed
   PT pt( nmr );
   tic = TIME_NOW();
   auto costPT = pt.solve( tmax, verbose );
   toc = ETS( tic );
   write_log( fid, "> timeoutPT ......... %d\n", pt.m_timeout ? 1 : 0 );
   write_log( fid, "> costPT ............ %d\n", costPT );
   write_log( fid, "> timePT (secs) ..... %ld\n", toc );

   fclose( fid );

   return EXIT_SUCCESS;
}