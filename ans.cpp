#pragma GCC optimize ("O3")
#pragma GCC target ("sse4")
#include <bits/stdc++.h>


//
#pragma comment(linker, "/stack:200000000")


bool ok(int x, int y) { return x >= 0 && y >= 0 && x < 3 && y < 3; }
int dx[4] = {1, 0, -1, 0}, dy[4] = {0, 1, 0, -1};

const int dx[8] = {1, 0, -1, 0, 1, 1, -1, -1}, dy[8] = {0, 1, 0, -1, -1, 1, -1, 1};

#define ll           long long
#define all(v)       v.begin(),v.end()
#define debug(...)   [](const auto&...x){ char c='='; cerr<<#__VA_ARGS__<<" "; ((cerr<<exchange(c,',')<<" "<<x),...); cerr<<endl; }(__VA_ARGS__);
#define nl           '\n'
#define sz(x)        (int)x.size()
#define ff           first
#define ss           second
#define convertToInt(x) stoi(x,nullptr,10)
#define vi           vector<int>


#define pb push_back
ll gcd(ll a, ll b) { return (b == 0LL ? a : gcd(b, a % b));}
ll lcm(ll a, ll b) { return (a * (b / gcd(a, b))); }
#define setbits(x)      __builtin_popcountll(x)
#define maxx(v)       *max_element(all(v));
#define minn(v)       *min_element(all(v));
#define rep(i,n)     for (int i = 0; i < n; i++)
#define repe(i,n)    for (int i = 0 ; i <= n ; i++)
#define rep1(i,n)    for (int i = 1 ; i < n ; i++)

#define rep1e(i,n)   for (int i = 1 ; i <= n ; i++)
#define repd(i,n)    for (int i = n - 1 ; i >= 0 ; i--)
#define repd1(i,n)   for (int i = n - 1 ; i >= 1 ; i--)
#define trav(val,v)  for (auto val : v)
#define travi(val,v) for (auto &val : v)
#define in(a)        cin >> a;
#define vi           vector<int>
#define vstr         vector<string>
#define vll          vector<ll>
#define sti          set<int>
#define stll         set<ll>
// mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

const int inf = 2e5 + 1;
const ll si = 1005;
const ll mod = 1e7 + 1;
const double pi = 3.141592653589793238;




bool ckmin(int& a, int b) { return b < a ? a = b, true : false; }

bool ckmax(int& a, int b) { return b > a ? a = b, true : false; }


// #ifdef LOCAL
// freopen("input.txt","r",stdin);
// freopen("output.txt","w",stdout);



#define to_check_adajcent_cells const int dx[4] = {1,0,-1,0}, dy[4] = {0,1,0,-1};
#define to_make_loop for(int i = 0; i < 4; i++){
#define index newX = x + dx[i]; newY = y + dy[i];
#define check_bounds bool ok(int x, int y) { return x >= 0 && y >= 0 && x < n && y < m; }

// the formula x | (1 << k) sets the kth bit of x to one,
//the formula x & ~(1 << k) sets the kth bit of x to zero, and
//the formula x ^ (1 << k) inverts the kth bit of x.
// The formula x & (x−1) sets the last one bit of x to zero, and
// the formula x & −x sets all the one bits to zero, except for the last one bit.
//  The formula x | (x−1) inverts all the bits after the last one bit.
// Also note that a positive number x is a
// power of two exactly when x & (x−1) = 0.


// few useful funtions rotate,bitset,distance,transform,find,count
// transform(a.begin(),a.end(),a.begin(),::toupper);
// is_sorted(begin,end) : used to check if arr is sorted bw begin & end

/*
 __builtin_clz(x):  the number of zeros at the beginning of the bit representation
 __builtin_ctz(x):  the number of zeros at the end of the bit representation
 __builtin_popcount(x):  the number of ones in the bit representation
 __builtin_parity(x):  the parity (even or odd) of the number of ones in the
bit representation
*/


void cout_the_subsequence(int arr[], int n) {
    int size = pow(2, n);
    for (int counter = 1; counter < size; counter++) {
        for (int i = 0; i < n; i++) {
            if (counter & (1 << i)) {
                cout << arr[i] << " ";
            }
        }
        cout << endl;
    }

}


// In ascii A = 65; Z = 90; a  = 97 ; z  = 122 ;  0 = 48  and 9 = 57;
// There's also __builtin_ffs(x) (Find  First Set) which returns (the index of the least significant bits of x) + 1.
// There's also __lg(x) which returns the index of the highest set bit.


ll gcd(ll a, ll b) { return (b == 0LL ? a : gcd(b, a % b));}
ll lcm(ll a, ll b) { return (a * (b / gcd(a, b))); }
const int dx[8] = {1, 0, -1, 0, 1, 1, -1, -1};
const int dy[8] = {0, 1, 0, -1, -1, 1, -1, 1};

bool  ok(int x, int y) { return x >= 0 && y >= 0 && x < n && y < m; }


class Graph {
public:
    map<int, bool> visited;
    map<int, list<int> > adj;

    // function to add an edge to graph
    void addEdge(int v, int w);
    queue<int> q;


    // DFS traversal of the vertices
    // reachable from v
    void printgraph(int v);
    void DFS(int v);
    void BFS(int v);
};

void Graph::addEdge(int u, int w) {
    adj[u].push_back(w);
}


void Graph::DFS(int v) {
    //Mark present node as visited
    visited[v] = true;
    cout << v << " ";
    for (auto edge : adj[v]) {
        if (!visited[edge]) {
            DFS(edge);
        }
    }
}

void Graph::BFS(int v) {
    cout << v << " ";
    visited[v] = true;
    for (auto val : adj[v]) {
        q.push(val);
        int temp = q.front();
        BFS(temp);
        q.pop();
    }
}


void Graph::printgraph(int V) {
    for (int i = 0 ; i < V ; i++) {
        for (auto val : adj[i]) {
            cout << val << " ";
        }
        cout << endl;
    }
}


function<int(int, int)> gcd = [&](int a, int b) {
    return b == 0 ? a : gcd(b , a % b);
};


template <typename T>
std::ostream& operator<<(std::ostream& out, const std::vector<T>& vec) {
    if (vec.empty()) {
        out << "[]";
        return out;
    }
    out << '[';
    for (int i = 0; i < vec.size() - 1; i++) {
        out << vec[i] << ", ";
    }
    return out << vec.back() << ']';
}

template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& out, const std::pair<T1, T2>& pair) {
    return out << '(' << pair.first << ", " << pair.second << ')';
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::deque<T>& deq) {
    if (deq.empty()) {
        out << "[]";
        return out;
    }
    out << '[';
    for (int i = 0; i < deq.size() - 1; i++) {
        out << deq[i] << ", ";
    }
    return out << deq.back() << ']';
}

template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& out, const std::unordered_map<T1, T2>& map) {
    out << '{';
    for (auto it = map.begin(); it != map.end(); it++) {
        std::pair<T1, T2> element = *it;
        out << element.first << ": " << element.second;
        if (std::next(it) != map.end()) {
            out << ", ";
        }
    }
    return out << '}';
}

template <typename T1, typename T2>
std::ostream& operator<<(std::ostream& out, const std::map<T1, T2>& map) {
    out << '{';
    for (auto it = map.begin(); it != map.end(); it++) {
        std::pair<T1, T2> element = *it;
        out << element.first << ": " << element.second;
        if (std::next(it) != map.end()) {
            out << ", ";
        }
    }
    return out << '}';
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::unordered_set<T>& set) {
    out << '{';
    for (auto it = set.begin(); it != set.end(); it++) {
        T element = *it;
        out << element;
        if (std::next(it) != set.end()) {
            out << ", ";
        }
    }
    return out << '}';
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::multiset<T>& set) {
    out << '{';
    for (auto it = set.begin(); it != set.end(); it++) {
        T element = *it;
        out << element;
        if (std::next(it) != set.end()) {
            out << ", ";
        }
    }
    return out << '}';
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::unordered_multiset<T>& set) {
    out << '{';
    for (auto it = set.begin(); it != set.end(); it++) {
        T element = *it;
        out << element;
        if (std::next(it) != set.end()) {
            out << ", ";
        }
    }
    return out << '}';
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const std::set<T>& set) {
    out << '{';
    for (auto it = set.begin(); it != set.end(); it++) {
        T element = *it;
        out << element;
        if (std::next(it) != set.end()) {
            out << ", ";
        }
    }
    return out << '}';
}

// Source: https://stackoverflow.com/a/31116392/12128483
template<typename Type, unsigned N, unsigned Last>
struct TuplePrinter {
    static void print(std::ostream& out, const Type& value) {
        out << std::get<N>(value) << ", ";
        TuplePrinter < Type, N + 1, Last >::print(out, value);
    }
};

template<typename Type, unsigned N>
struct TuplePrinter<Type, N, N> {
    static void print(std::ostream& out, const Type& value) {
        out << std::get<N>(value);
    }
};

template<typename... Types>
std::ostream& operator<<(std::ostream& out, const std::tuple<Types...>& value) {
    out << '(';
    TuplePrinter < std::tuple<Types...>, 0, sizeof...(Types) - 1 >::print(out, value);
    return out << ')';
}


/*
#Sigma Rule :-
        If you are not finding any mistake in code then there is mistake in your logic.

*/

// Too dumb to exit
/*
////////////////////////////////////////////////
///  Rishabh Kumar Pandey                    ///
///  Vellore Institute Of Technology,Vellore ///
////////////////////////////////////////////////
*/


// freopen("piggyback.in", "r", stdin);
// the following line creates/overwrites the output file
// freopen("piggyback.out", "w", stdout);


vector<bool> isprime(201, true);
void pre() {
    isprime[0] = false;
    isprime[1] = false;
    for (int i = 2 ; i <= 201 ; i++) {
        if (isprime[i]) {
            for (int j = i * i ; j <= 201 ; j += i ) {
                isprime[j] = false;
            }
        }
    }
}



407833339917496991 742004829320395233


334171489402898242



// for the implementation of lambda functions which involve recursion
namespace std {

template<class Fun>
class y_combinator_result {
    Fun fun_;
public:
    template<class T>
    explicit y_combinator_result(T &&fun): fun_(std::forward<T>(fun)) {}

    template<class ...Args>
    decltype(auto) operator()(Args &&...args) {
        return fun_(std::ref(*this), std::forward<Args>(args)...);
    }
};

template<class Fun>
decltype(auto) y_combinator(Fun &&fun) {
    return y_combinator_result<std::decay_t<Fun>>(std::forward<Fun>(fun));
}

}


// LIS in nlong(n)
int lis(vector<int> const& a) {
    int n = a.size();
    const int INF = 1e9;
    vector<int> d(n+1, INF);
    d[0] = -INF;

    for (int i = 0; i < n; i++) {
        int l = upper_bound(d.begin(), d.end(), a[i]) - d.begin();
        if (d[l-1] < a[i] && a[i] < d[l])
            d[l] = a[i];
    }

    int ans = 0;
    for (int l = 0; l <= n; l++) {
        if (d[l] < INF)
            ans = l;
    }
    return ans;
}
