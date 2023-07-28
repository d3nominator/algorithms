// File to be included
// Common file
#include <ext/pb_ds/assoc_container.hpp>
// Including tree_order_statistics_node_update
#include <ext/pb_ds/tree_policy.hpp>
#include <ext/pb_ds/detail/standard_policies.hpp>

//
// visit this link to study it in a proper way :
// https://codeforces.com/blog/entry/11080


using namespace __gnu_pbds;
#define ordered_set  tree<int, null_type, less_equal<int>, rb_tree_tag, tree_order_statistics_node_update> 
// use less_equal for multiset and less for set ^


// Find by order and order of key
// cout<< *X.find_by_order(k)    <<endl;  //  K-th element in a set (counting from zero)
// cout<< X.order_of_key(k)   <<endl; //  Number of items strictly smaller than k.


William lin version
template<class T> using oset =  tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>; 
oset<array<int,2>> s;


#include <ext/pb_ds/assoc_container.hpp>
// Including tree_order_statistics_node_update
#include <ext/pb_ds/tree_policy.hpp>
#include <ext/pb_ds/detail/standard_policies.hpp>


using namespace __gnu_pbds;
#define ordered_set  tree<int, null_type, less_equal<int>, rb_tree_tag, tree_order_statistics_node_update> 



// Include the below one ->
#include <ext/pb_ds/assoc_container.hpp>
// Including tree_order_statistics_node_update
#include <ext/pb_ds/tree_policy.hpp>
#include <ext/pb_ds/detail/standard_policies.hpp>

using namespace __gnu_pbds;
#define ordered_set  tree<int, null_type, less_equal<int>, rb_tree_tag, tree_order_statistics_node_update> 















///////////////////////////////////////////Ordered multiset implementation
// File to be included
// Common file
#include <ext/pb_ds/assoc_container.hpp>
// Including tree_order_statistics_node_update
#include <ext/pb_ds/tree_policy.hpp>
#include <ext/pb_ds/detail/standard_policies.hpp>
using namespace __gnu_pbds;
using ll = long long;
struct ordered_multiset
{ // multiset supporting duplicating values in set
    ll len = 0;
    const ll ADD = 1000010;
    const ll MAXVAL = 1000000010;
    map<ll, ll> mp; // hash = 96814
    tree<ll, null_type, less<ll>, rb_tree_tag, tree_order_statistics_node_update> T;

    ordered_multiset()
    {
        len = 0;
        T.clear(), mp.clear();
    }

    inline void insert(ll x)
    {
        len++, x += MAXVAL;
        ll c = mp[x]++;
        T.insert((x * ADD) + c);
    }

    inline void erase(ll x)
    {
        x += MAXVAL;
        ll c = mp[x];
        if (c)
        {
            c--, mp[x]--, len--;
            T.erase((x * ADD) + c);
        }
    }

    inline ll kth(ll k)
    { // 1-based index,  returns the
        if (k < 1 || k > len)
            return -1;                  // K'th element in the treap,
        auto it = T.find_by_order(--k); // -1 if none exists
        return ((*it) / ADD) - MAXVAL;
    }

    inline ll lower_bound(ll x)
    { // Count of value < x in treap
        x += MAXVAL;
        ll c = mp[x];
        return (T.order_of_key((x * ADD) + c));
    }

    inline ll upper_bound(ll x)
    { // Count of value <= x in treap
        x += MAXVAL;
        ll c = mp[x];
        return (T.order_of_key((x * ADD) + c));
    }

    inline ll size() { return len; } // Number of elements in treap
};