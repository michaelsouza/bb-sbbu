#include <algorithm>
#include <chrono>
#include <climits>
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
      return ( m_i + 3 <= s.m_i ) && ( s.m_j <= m_j );
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

weight_t order_sbbu( NMR& nmr, std::vector<int>& order ) {
   auto& E = nmr.m_E;

   order.clear();
   for ( auto&& kv : E )
      order.push_back( kv.first ); // list of edges eid

   std::sort( order.begin(), order.end(),
       [ &E ]( const int& eidA, const int& eidB ) -> bool {
         const auto& eA = E[ eidA ];
         const auto& eB = E[ eidB ];
          if ( eA.m_j == eB.m_j )
             return eA.m_i < eB.m_i;
          return eA.m_j < eB.m_j;
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

class BBPerm {
private:
   void setMembers() {
      m_state = 'n';
      m_idx = -1;
      m_order.resize( m_keys.size() );
      for ( auto&& key : m_keys )
         m_bst.add( key );
   }

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
      setMembers();
   }

   BBPerm( std::vector<int>& keys ) {
      std::copy( keys.begin(), keys.end(), m_keys.begin() );
      setMembers();
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
            auto key = m_bst.minKey( true );
            m_order[ ++m_idx ] = key;
            return key;
         }
         m_state = 'p';
         return next();
      }

      if ( m_idx == -1 )
         return -1;

      const auto keyOld = m_order[ m_idx ];
      auto key = m_bst.minKeyGT( keyOld, true );
      m_bst.add( keyOld );
      m_order[ m_idx ] = key;
      if ( key < 0 ) {
         --m_idx;
         return next();
      }
      m_state = 'n';
      return key;
   }

   void prune() {
      m_state = 'p';
   }
};

class BB {
public:
   int m_idx;
   NMR& m_nmr;
   BBPerm m_perm;
   bool m_timeout;
   size_t m_nedges;
   std::vector<int> m_order;
   std::map<int, NMREdge>& m_E;
   std::map<int, NMRSegment>& m_S;
   // array of the cost added by each eid
   std::vector<weight_t> m_c;

   BB( NMR& nmr )
       : m_nmr( nmr ), m_E( nmr.m_E ), m_S( nmr.m_S ), m_perm( nmr.m_pruneEdges ) {
      m_idx = -1;
      m_nedges = m_E.size();
      m_order.resize( m_nedges );
      m_timeout = false;
      m_c.resize( m_nedges );
      std::fill( m_c.begin(), m_c.end(), 0 );
   }

   /**
    * @brief
    *
    * @param C C[id] is the number of edges on m_order covering the segment sid.
    * @param U set of available segments
    * @return weight_t
    */
   weight_t order_rem( std::vector<int>& C, std::set<int>& U, weight_t& costRLX ) {
      weight_t total_cost = 0;
      while ( m_idx >= m_perm.m_idx && m_perm.m_order[ m_idx ] != m_order[ m_idx ] ) {
         auto eid = m_order[ m_idx ];
         // set invalid value
         m_order[ m_idx ] = -1;
         total_cost += m_c[ m_idx ];
         const auto& e = m_E[ eid ];
         for ( auto&& sid : e.m_sid ) {
            --C[ sid ];
            if ( C[ sid ] == 0 ) {
               U.insert( sid );
               const auto s = m_S[ sid ];
               // increment the relaxed cost
               costRLX += s.m_weight;
            }
         }
         --m_idx;
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
   weight_t order_add( int eid, std::vector<int>& C, std::set<int>& U, weight_t& costRLX ) {
      m_order[ ++m_idx ] = eid;
      weight_t eid_cost = 1;
      const auto& e = m_E[ eid ];
      for ( auto&& sid : e.m_sid ) {
         ++C[ sid ];
         // the current eid is the only one covering the sid
         if ( C[ sid ] == 1 ) {
            const auto& s = m_S[ sid ];
            eid_cost *= s.m_weight;
            U.erase( sid );
            // reduce the relaxed cost
            costRLX -= s.m_weight;
         }
      }
      eid_cost = eid_cost > 1 ? eid_cost : 0;
      m_c[ m_idx ] = eid_cost;
      return eid_cost;
   }

   weight_t solve( size_t tmax = 3600 ) {
      auto tic = TIME_NOW();
      weight_t costUB = WEIGHT_MAX;
      std::fill( m_order.begin(), m_order.end(), -1 );
      std::vector<int> orderOPT;
      costUB = order_sbbu( m_nmr, orderOPT );

      // C[sid] : number of edges already included in the order that cover segment sid
      std::vector<int> C( m_S.size() + 1 ); // the segments are one-based
      std::fill( C.begin(), C.end(), 0 );

      // U: set of the uncovered segments
      std::set<int> U;
      for ( auto&& kv : m_S )
         U.insert( kv.first );

      // first cost_relax
      auto costRELAX = cost_relax( U, m_S );

      // solution found
      if ( costRELAX == costUB ) {
         std::copy( orderOPT.begin(), orderOPT.end(), m_order.begin() );
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

      int eid = m_perm.next();
      // loop through all permutations
      while ( eid > -1 ) {
         auto toc = ETS( tic );
         if ( toc > tmax ) {
            m_timeout = true;
            printf( "> timeoutBB %ld seconds\n", toc );
            break;
         }

         costACC -= order_rem( C, U, costRLX );
         costEID = order_add( eid, C, U, costRLX );

         costACC += costEID;
         costLB = costACC + costRLX;
         
         // sanity check (something went wrong)
         if( costLB < costRELAX ){
            costUB = costLB;
            std::copy( m_order.begin(), m_order.end(), orderOPT.begin() );
            break;
         }

         // update solution
         if ( ( costRLX == 0 ) && ( costLB < costUB ) ) {
            costUB = costLB;
            std::copy( m_order.begin(), m_order.end(), orderOPT.begin() );
            if ( costRELAX >= costLB ) break;
         }

         // prune when
         // 1) costLB > costUB: the prefix will not generate a improved solution
         // 2) costRLX == 0: there is no remaining segment to be solved
         // 3) costEID == 0: push the solved eid to the end of the permutation
         if ( ( costLB > costUB ) || ( costRLX == 0 ) || ( costEID == 0 ) ) m_perm.prune();

         // TODO Jump permutation when U is empty
         // TODO Consider only the eid's covering the sid's on U
         eid = m_perm.next();
      }

      // sanity check (something went wrong)
      if ( costRELAX > costUB )
         throw std::runtime_error( "Unfeasible cost found (costRELAX < costUB)." );

      std::copy( orderOPT.begin(), orderOPT.end(), m_order.begin() );
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

   fflush(fid);
}

int call_solvers( int argc, char* argv[] ) {
   std::string fnmr = "/home/michael/gitrepos/bb-sbbu/DATA_TEST/testC.nmr";
   size_t tmax = 3600;
   bool clean_log = false;
   for ( int i = 0; i < argc; ++i ) {
      auto arg = argv[ i ];
      if ( strcmp( arg, "-fnmr" ) == 0 ) fnmr = argv[ i + 1 ];
      if ( strcmp( arg, "-tmax" ) == 0 ) tmax = atoi( argv[ i + 1 ] );
      if ( strcmp( arg, "-clean_log" ) == 0 ) clean_log = true;
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

   write_log( fid, "> clean_log ......... %d\n", clean_log ? 1 : 0 );
   write_log( fid, "> tmax (secs) ....... %ld\n", tmax );
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
   write_log( fid, "> timeSBBU (secs) ... %ld\n", toc );

   // call order_bb if needed
   BB bb( nmr );
   tic = TIME_NOW();
   auto costBB = bb.solve( tmax = tmax );
   toc = ETS( tic );
   write_log( fid, "> timeoutBB ......... %d\n", bb.m_timeout ? 1 : 0 );
   write_log( fid, "> costBB ............ %d\n", costBB );
   write_log( fid, "> timeBB (secs) ..... %ld\n", toc );

   fclose( fid );

   return EXIT_SUCCESS;
}