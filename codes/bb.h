#include <algorithm>
#include <chrono>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <regex>
#include <set>
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
   std::set<int> m_eid;

   NMRSegment() {
      m_sid = -1;
      m_i = -1;
      m_j = -1;
      m_weight = -1;
   }

   NMRSegment( const NMRSegment& s )
       : m_eid( s.m_eid ) {
      m_sid = s.m_sid;
      m_i = s.m_i;
      m_j = s.m_j;
      m_weight = s.m_weight;
   }

   NMRSegment( int i, int j ) {
      m_sid = ++NMRSegment::SID;
      m_i = i;
      m_j = j;
      update_weight();
   }

   static void resetSID() {
      NMRSegment::SID = 0;
   }

   void add_eid( int eid ) {
      m_eid.insert( eid );
   }

   void update_weight() {
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
   std::set<int> m_sid;

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

   void add_sid( int sid ) {
      m_sid.insert( sid );
   }

   bool check_cover( NMRSegment& s ) {
      return (m_i + 3 <= s.m_i) && (s.m_j <= m_j);
   }

   static void resetEID() {
      NMREdge::EID = 0;
   }
};

int NMREdge::EID;

class NMR {
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
      for ( auto&& i : I ) E[i] = {};

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
            if ( e.check_cover( s ) ) {
               e.add_sid( s.m_sid );
               s.add_eid( e.m_eid );
            }
   }

   void setOrderingData() {
      for ( auto&& e : m_pruneEdges )
         m_E[ e.m_eid ] = e;
      for ( auto&& s : m_segments )
         m_S[ s.m_sid ] = s;
   }

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
      
      // open fnmr file
      std::ifstream fid( fnmr.c_str() );
      if(!fid.good()) throw std::runtime_error("Could not open fnmr file.");

      int i, j;
      NMREdge::resetEID(); // set EID=0
      while ( fid >> i >> j ) m_edges.push_back( NMREdge( i, j ) );
      
      m_nnodes = 0;
      for ( auto&& e : m_edges ) {
         if ( m_nnodes < e.m_j ) m_nnodes = e.m_j;
         if ( e.m_j > e.m_i + 3 ) m_pruneEdges.push_back( e );
      }
      setSegments();
      setOrderingData();
   }
};

weight_t order_cost( std::vector<int>& order, std::map<int, NMREdge>& E, std::map<int, NMRSegment>& S, weight_t costUB = WEIGHT_MAX ) {
   weight_t total_cost = 0; // total cost
   // sid is added to B by the first edge that covers it.
   std::set<int> B;   
   for ( auto&& eid : order ) {
      weight_t edge_cost = 1;
      const auto& e = E[ eid ];
      for ( auto&& sid : e.m_sid ) {
         // sid not in S, it does not need to be covered
         if( S.find( sid ) == S.end() ) continue;
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

weight_t order_sbbu( NMR& nmr, std::vector<int>& order ) {
   auto& E = nmr.m_E;

   order.clear();
   for ( auto&& kv : E )
      order.push_back( kv.first ); // list of edges eid

   std::sort( order.begin(), order.end(),
       [ &E ]( const int& eidA, const int& eidB ) -> bool {
          if ( E[ eidA ].m_j == E[ eidB ].m_j )
             return E[ eidA ].m_i < E[ eidB ].m_j;
          return E[ eidA ].m_j < E[ eidB ].m_j;
       } );

   return order_cost( order, E, nmr.m_S );
}

weight_t order_brute( NMR& nmr, std::vector<int>& orderOPT ) {
   auto& E = nmr.m_E;
   auto& S = nmr.m_S;
   weight_t costOPT = WEIGHT_MAX;

   std::vector<int> order;
   for ( auto&& kv : E )
      order.push_back( kv.first );

   orderOPT.resize( order.size() );

   do {
      auto c = order_cost( order, E, S );
      if ( c < costOPT ) {
         costOPT = c;
         std::copy( order.begin(), order.end(), orderOPT.begin() );
      }
   } while ( std::next_permutation( order.begin(), order.end() ) );

   return costOPT;
}

weight_t cost_relax( std::set<int>& U, std::map<int, NMRSegment>& S ) {
   weight_t total_cost = 0;
   for ( auto&& sid : U ) total_cost += S[ sid ].m_weight;
   return total_cost;
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

   BSTNode* min( BSTNode* x = nullptr ) {
      if ( x == nullptr ) x = m_root;
      while ( x->m_lft != nullptr ) x = x->m_lft;
      return x;
   }

   BSTNode* max() {
      BSTNode* x = m_root;
      while ( x->m_rht != nullptr ) x = x->m_rht;
      return x;
   }

   BSTNode* minGT( int key ) {
      BSTNode* u = nullptr;
      BSTNode* x = m_root;
      int minKey = INT_MAX;
      while ( x != nullptr ) {
         if ( x->m_key > key ) {
            minKey = key;
            u = x;
            x = x->m_lft;
         }
         else { // x->m_key >= key
            x = x->m_rht;
         }
      }
      return u;
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

   void rem( BSTNode* z, bool delz = false ) {
      --m_size;
      if ( z->m_lft == nullptr ) { // (possibly) only rht child
         transplant( z, z->m_rht );
      }
      else if ( z->m_rht == nullptr ) { // only lft child
         transplant( z, z->m_lft );
      }
      else {
         auto y = min();
         if ( y != z->m_rht ) {
            transplant( y, y->m_rht );
            y->m_rht = z->m_rht;
            y->m_rht->m_p = y;
         }
         transplant( z, y );
         y->m_lft = z->m_lft;
         y->m_lft->m_p = y;
      }
      if ( delz )
         delete z;
   }
};

class BBPerm {
public:
   std::vector<int> m_keys;
   // current order
   std::vector<int> m_order;
   // index of the element last element inserted
   int m_idx;
   // n:normal, p:prune
   char m_state;
   BST m_bst;

   BBPerm( std::vector<NMREdge>& E ) {
      for ( auto&& e : E )
         m_keys.push_back( e.m_eid );
      m_state = 'n';
      m_idx = -1;
      m_order.resize( m_keys.size() );
   }

   BBPerm( std::vector<int>& keys ) {
      std::copy( keys.begin(), keys.end(), m_keys.begin() );
      m_state = 'n';
      m_idx = -1;
      m_order.resize( m_keys.size() );
   }

   /**
    * @brief Return the index of the last component in the current order. If there is no
    *        next component, the index -1 is returned.
    *
    * @return int
    */
   int next() {
      if ( m_state == 'n' ) {
         if ( m_bst.m_size > 0 ) {
            BSTNode* emin = m_bst.min();
            int key = emin->m_key;
            m_bst.rem( emin );
            m_order[ ++m_idx ] = emin->m_key;
            return emin->m_key;
         }
         m_state = 'p';
         return next();
      }

      if ( m_idx == -1 )
         return -1;

      int key = m_order[ m_idx ];
      BSTNode* emin = m_bst.minGT( key );
      m_bst.add( key );
      if ( emin == nullptr ) {
         m_order[ m_idx ] = -1;
         --m_idx;
         return next();
      }
      m_order[ m_idx ] = emin->m_key;
      m_state = 'n';
      return key;
   }

   void prune() {
      m_state = 'p';
   }
};

class BB {
public:
   NMR& m_nmr;
   std::map<int, NMREdge>& m_E;
   std::map<int, NMRSegment>& m_S;
   size_t m_nedges;
   int m_idx;
   BBPerm m_perm;
   std::vector<int> m_order;
   bool m_timeout;

   BB( NMR& nmr )
       : m_nmr( nmr ), m_E( nmr.m_E ), m_S( nmr.m_S ), m_perm( nmr.m_pruneEdges ) {
      m_idx = -1;
      m_nedges = m_E.size();
      m_order.resize( m_nedges );
      m_timeout = false;
   }

   /**
    * @brief
    *
    * @param C C[id] is the number of edges on m_order covering the segment sid.
    * @param U set of available segments
    * @return weight_t
    */
   weight_t order_rem( std::map<int, int>& C, std::set<int>& U ) {
      auto idx = m_idx;
      weight_t total_cost = 0;
      while ( idx >= m_perm.m_idx && m_perm.m_order[ idx ] != m_order[ idx ] ) {
         auto eid = m_order[ idx ];
         // set invalid value
         m_order[ idx ] = -1;
         weight_t eid_cost = 1;
         for ( auto&& sid : m_E[ eid ].m_sid ) {
            --C[ sid ];
            if ( C[ sid ] == 0 ) {
               eid_cost *= m_S[ sid ].m_weight;
               U.insert( sid );
            }
         }
         if ( eid_cost > 1 )
            total_cost += eid_cost;
         --idx;
      }
      return total_cost;
   }

   /**
    * @brief Add edge eid to the m_order and update C and U.
    *
    * @param eid
    * @param C C[sid] : number of edges on m_order covering segment sid
    * @param U set of the uncovered segments.
    * @return weight_t
    */
   weight_t order_add( int eid, std::map<int, int>& C, std::set<int>& U ) {
      m_order[ m_perm.m_idx ] = eid;
      weight_t eid_cost = 1;
      for ( auto&& sid : m_E[ eid ].m_sid ) {
         ++C[ sid ];
         // the current eid is the only one covering the sid
         if ( C[ sid ] == 1 ) {
            eid_cost *= m_S[ sid ].m_weight;
            U.erase( sid );
         }
      }
      return eid_cost > 1 ? eid_cost : 0;
   }

   weight_t solve( bool unpickling = false, float tmax = 60 ) {
      auto tic = TIME_NOW();
      weight_t costUB = WEIGHT_MAX;
      std::vector<int> orderOPT;
      if ( unpickling )
         throw std::invalid_argument( "not implemented." );
      else // initial optimal solution
         costUB = order_sbbu( m_nmr, m_order );

      // C[sid] : number of edges already included in the order that cover segment sid
      std::map<int, int> C;
      for ( auto&& kv : m_S ) C[ kv.first ] = 0;

      // U: set of the uncovered segments
      std::set<int> U;
      for ( auto&& kv : C ) U.insert( kv.first );

      // first cost_relax
      auto costLB = cost_relax( U, m_S );
      if ( costLB == costUB ) return costUB;

      weight_t partial_cost = 0;
      int eid = m_perm.next();
      // loop through all permutations
      while ( eid > -1 ) {
         partial_cost -= order_rem( C, U );
         partial_cost += order_add( eid, C, U );
         m_idx = m_perm.m_idx;
         // when U is empty, the partial_cost is total.
         costLB = partial_cost + cost_relax( U, m_S );
         auto toc = ETS( tic );
         if ( toc > tmax ) {
            m_timeout = true;
            printf( "> timeoutBB %ld seconds\n", toc );
            // m_dump(C, U)
            break;
         }
         if ( costLB >= costUB )
            m_perm.prune();
         else if ( m_perm.m_idx == ( m_nedges - 1 ) && costLB < costUB ) {
            costUB = costLB;
            std::copy( m_order.begin(), m_order.end(), orderOPT.begin() );
         }
         eid = m_perm.next();
      }
      std::copy( orderOPT.begin(), orderOPT.end(), m_order.begin() );
      return costUB;
   }
};

bool exists( std::string fname )
{
   std::ifstream fid( fname.c_str() );
   return fid.good();
}

void write_log( FILE* fid, const char* fmt, ... )
{
   va_list args;
   va_start( args, fmt );
   vprintf( fmt, args );
   vfprintf( fid, fmt, args );
   va_end( args );
}

int call_solvers( int argc, char* argv[] )
{
   std::string fnmr = "/home/michael/gitrepos/bb-sbbu/DATA_TEST/testC.nmr";
   int tmax = 1;
   bool clean_log = false;
   for ( int i = 0; i < argc; ++i )
   {
      auto arg = argv[ i ];
      if ( arg == "-fnmr" ) fnmr = argv[ i + 1 ];
      if ( arg == "-tmax" ) tmax = atof( argv[ i + 1 ] );
      if ( arg == "-clean_log" ) clean_log = true;
   }

   auto flog = std::regex_replace( fnmr, std::regex( ".nmr" ), ".log" );
   // check if already has a log file
   if ( !clean_log && exists( flog ) )
   {
      printf( "> skip (already solved) %s\n", fnmr.c_str() );
      return EXIT_SUCCESS;
   }
   // create log file
   FILE* fid = fopen( flog.c_str(), "w" );
   write_log( fid, "> fnmr %s\n", fnmr.c_str() );

   // read instance
   NMR nmr( fnmr );
   auto& E = nmr.m_E;
   auto& S = nmr.m_S;

   write_log( fid, "> tmax (secs) ....... %g\n", tmax );
   write_log( fid, "> nnodes ............ %d\n", nmr.m_nnodes );
   write_log( fid, "> lenE .............. %d\n", E.size() );
   write_log( fid, "> lenS .............. %d\n", S.size() );

   std::set<int> U;
   for ( auto&& kv : S ) U.insert( kv.first );

   auto costRELAX = cost_relax( U, S );
   write_log( fid, "> costRELAX ......... %d\n", costRELAX );

   // call order_sbbu
   std::vector<int> orderSBBU;
   auto tic = TIME_NOW();
   auto costSBBU = order_sbbu( nmr, orderSBBU );
   auto toc = ETS( tic );
   write_log( fid, "> costSBBU .......... %d\n", costSBBU );
   write_log( fid, "> timeSBBU (secs) ... %g\n", toc );

   // call order_bb if needed
   BB bb( nmr );
   tic = TIME_NOW();
   auto costBB = bb.solve( tmax = tmax );
   toc = ETS( tic );
   write_log( fid, "> timeoutBB ......... %s\n", bb.m_timeout ? "true" : "false" );
   write_log( fid, "> costBB ............ %d\n", costBB );
   write_log( fid, "> timeBB (secs) ..... %g\n", toc );

   fclose( fid );

   return EXIT_SUCCESS;
}