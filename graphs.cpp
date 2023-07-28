// DFS implementation
adjacency list representation : vector<int> adj[N];
bool visited[N];
void dfs(int s) {
    if (visited[s]) return;
    visited[s] = true;
    // process node s
    for (auto u : adj[s]) {
        dfs(u);
    }
}

// graph class for easy implemetation
class Graph {
    int V;
    list<int>* adj;
public:
    Graph(int V)
    {
        this->V = V;
        adj = new list<int>[V];
    }
    void addEdge(int v, int w)
    {
        adj[v].push_back(w);
        adj[w].push_back(v);
    }
    void DFSUtil(int v, bool visited[], vi &v1)
    {
        visited[v] = true;
        v1.pb(v);
        list<int>::iterator i;
        for (i = adj[v].begin(); i != adj[v].end(); ++i)
            if (!visited[*i])
                DFSUtil(*i, visited, v1);
    }
    void connectedComponents(vector<vi> &p)
    {
        bool* visited = new bool[V];
        for (int v = 0; v < V; v++)
            visited[v] = false;
        for (int v = 0; v < V; v++) {
            if (visited[v] == false) {
                vi v1;
                DFSUtil(v, visited, v1);
                total++;
                p.push_back(v1);
            }
        }
    }
};



//struct for DSU in graphs
struct dsu
{
    vector<int> leader, siz, edg;
    void init(int n)
    {
        leader.resize(n + 1);
        siz.resize(n + 1, 1);
        edg.resize(n + 1, 0);
        for (int i = 1; i <= n; i++) leader[i] = i;
    }
    int get(int curr)
    {
        if (leader[curr] == curr) return curr;
        return leader[curr] = get(leader[curr]);
    }
    void merg(int a, int b)
    {
        if (get(a) != get(b)) {
            siz[get(b)] += siz[get(a)];
            edg[get(b)] += edg[get(a)] + 1;
            leader[get(a)] = get(b);
        }
        else edg[get(b)]++;
    }
};



// UFDS
struct UnionFind {
    //If the vertex is a parent, this value will be the number of vertices that belongs to the set, multiplied by -1
    //Otherwise, this value will be the parent's id
    vector<int> r;
    UnionFind(int N) {
        r = vector<int>(N, -1);
    }
    int root(int x) {
        if (r[x] < 0) return x;
        return r[x] = root(r[x]);
    }
    bool unite(int x, int y) {
        x = root(x);
        y = root(y);
        if (x == y) return false;
        if (r[x] > r[y]) swap(x, y);
        r[x] += r[y];
        r[y] = x;
        return true;
    }
    int size(int x) {
        return -r[root(x)];
    }
};



// Union find for finding the size of sets
// use the following
int parent[5000001];
int siz[5000001];

void make_set(int v) {
    parent[v] = v;
    siz[v] = 1LL;
}

int find_set(int v) {
    if (v == parent[v])
        return v;
    return parent[v] = find_set(parent[v]);
}

void union_sets(int a, int b) {
    a = find_set(a);
    b = find_set(b);
    if (a != b) {
        if (siz[a] < siz[b])
            swap(a, b);
        parent[b] = a;
        siz[a] += siz[b];
    }
}

++++ +


for (int i = 1LL ; i <= n ; i++) {
    make_set(i);
}




// Jiangly version of dsu
struct DSU {
    std::vector<int> f, siz;
    DSU(int n) : f(n), siz(n, 1) { std::iota(f.begin(), f.end(), 0); }
    int leader(int x) {
        while (x != f[x]) x = f[x] = f[f[x]];
        return x;
    }
    bool same(int x, int y) { return leader(x) == leader(y); }
    bool merge(int x, int y) {
        x = leader(x);
        y = leader(y);
        if (x == y) return false;
        siz[x] += siz[y];
        f[y] = x;
        return true;
    }
    int size(int x) { return siz[leader(x)]; }
};


//Segmentree
struct Node {
    Node() {
    }
    Node(char c) {
        ca = cab = cbc = cabc = pa = pb = pc = 0;
        if (c == 'a') ca++, pa++;
        if (c == 'b') pb++;
        if (c == 'c') pc++;
    }
    ll ca, cab, cbc, cabc, pa, pb, pc;
};

struct SegTree {
    int n;
    vector<int> l, r;
    vector<Node> a, tree;
    Node Null;
    SegTree(int size) {
        n = size;
        l.resize(4 * n);
        r.resize(4 * n);
        tree.resize(4 * n);
        a.resize(n + 1);
        Null = Node(0);
        init(1, n, 1);
    }

    void init(int lo, int hi, int node) {
        l[node] = lo;
        r[node] = hi;
        if (lo < hi) {
            init(lo, lo + (hi - lo) / 2, 2 * node);
            init(lo + (hi - lo) / 2 + 1, hi, 2 * node + 1);
        }
    }

    Node merge(Node &A, Node &B) {
        Node C;
        C.pa = A.pa + B.pa;
        C.pb = A.pb + B.pb;
        C.pc = A.pc + B.pc;
        C.ca = C.pa;
        C.cab = min(A.pa + B.cab, B.pb + A.cab);
        C.cbc = min(A.pb + B.cbc, B.pc + A.cbc);
        C.cabc = min(A.cab + B.cbc, min(A.pa + B.cabc, B.pc + A.cabc));
        return C;
    }

    void create(int node) {
        if (l[node] == r[node])
            tree[node] = a[l[node]];
        else {
            create(2 * node);
            create(2 * node + 1);
            tree[node] = merge(tree[2 * node], tree[2 * node + 1]);
        }
    }

    Node query(int lo, int hi, int node) {
        if (l[node] > hi || r[node] < lo) return Null;
        else if (l[node] >= lo && r[node] <= hi) return tree[node];
        else {
            Node A = query(lo, hi, 2 * node);
            Node B = query(lo, hi, 2 * node + 1);
            Node ans = merge(A, B);
            return ans;
        }

    }

    void update(int pos, int node) {
        if (l[node] <= pos && r[node] >= pos) {
            if (l[node] == r[node]) tree[node] = a[pos];
            else {
                update(pos, 2 * node);
                update(pos, 2 * node + 1);
                tree[node] = merge(tree[2 * node], tree[2 * node + 1]);
            }
        }
        else return;
    }
};


// Top sort
template <bool oneindexed = true> struct topsort {
    vector<vector<int>> edges;
    vector<int> indeg;
    int n;

    topsort() {}

    topsort(int N) {
        n = N;
        edges.resize(n);
        indeg.resize(n);
    }

    void add_edge(int x, int y) {
        edges[x - oneindexed].push_back(y - oneindexed);
        indeg[y - oneindexed]++;
    }

    vector<int> order() {
        // Return a topological ordering, or an empty vector if there is none
        vector<int> res, cur;
        for (int i = 0; i < n; i++) {
            if (indeg[i] == 0) cur.push_back(i);
        }
        while (!cur.empty()) {
            int node = cur.back();
            cur.pop_back();
            res.push_back(node + oneindexed);
            for (int i : edges[node]) {
                if (--indeg[i] == 0) cur.push_back(i);
            }
        }
        if (res.size() < n) res.clear();
        return res;
    }

    vector<int> minorder() {
        // Return the lexicographically minimal topological ordering, or an empty vector if there is none
        vector<int> res;
        priority_queue<int, vector<int>, greater<int>> cur;
        for (int i = 0; i < n; i++) {
            if (indeg[i] == 0) cur.push(i);
        }
        while (!cur.empty()) {
            int node = cur.top();
            cur.pop();
            res.push_back(node + oneindexed);
            for (int i : edges[node]) {
                if (--indeg[i] == 0) cur.push(i);
            }
        }
        if (res.size() < n) res.clear();
        return res;
    }

    vector<int> maxorder() {
        // Return the lexicographically maximal topological ordering, or an empty vector if there is none
        vector<int> res;
        priority_queue<int> cur;
        for (int i = 0; i < n; i++) {
            if (indeg[i] == 0) cur.push(i);
        }
        while (!cur.empty()) {
            int node = cur.top();
            cur.pop();
            res.push_back(node + oneindexed);
            for (int i : edges[node]) {
                if (--indeg[i] == 0) cur.push(i);
            }
        }
        if (res.size() < n) res.clear();
        return res;
    }
};



//
struct UnionFind {
    int size;
    vector<int> parent;
    vector<int> rank;
    vector<ll> v, e;

    UnionFind() {}
    UnionFind(int size) {
        this->size = size;
        parent.resize(size + 1);
        rank.resize(size + 1);
        v.resize(size + 1);
        e.resize(size + 1);
        init();
    }
    void init() {
        for (int i = 0; i <= size; i++) {
            parent[i] = i, rank[i] = 0;
            v[i] = 1, e[i] = 0;
        }
    }
    int root(int i) {
        if (parent[i] == i) return i;
        return parent[i] = root(parent[i]);
    }
    bool same(int i, int j) {
        return root(i) == root(j);
    }
    void merge(int i, int j) { // j will become new root
        parent[i] = j;
        v[j] += v[i];
        e[j] += e[i] + 1;
    }
    void unite(int i, int j) {
        int root_i = root(i), root_j = root(j);
        if (root_i == root_j) {
            e[root_i]++;
            return;
        }
        if (rank[root_i] < rank[root_j]) merge(root_i, root_j);
        else merge(root_j, root_i);
        if (rank[root_i] == rank[root_j]) rank[root_i]++;
    }
};



//functional graphs

template <typename T = int>
struct FunctionalGraph {
    int N, M;
    vc<int> TO;
    vc<T> wt;
    vc<int> root;
    Graph<T, 1> G;

    FunctionalGraph() {}
    FunctionalGraph(int N) : N(N), M(0), TO(N, -1), wt(N), root(N, -1) {}

    void add(int a, int b, T c = 1) {
        assert(0 <= a && a < N);
        assert(TO[a] == -1);
        ++M;
        TO[a] = b;
        wt[a] = c;
    }

    Tree<Graph<T, 1>> build() {
        assert(N == M);
        UnionFind uf(N);
        FOR(v, N) if (!uf.merge(v, TO[v])) { root[v] = v; }
        FOR(v, N) if (root[v] == v) root[uf[v]] = v;
        FOR(v, N) root[v] = root[uf[v]];

        G.resize(N + 1);
        FOR(v, N) {
            if (root[v] == v)
                G.add(N, v, wt[v]);
            else
                G.add(TO[v], v, wt[v]);
        }
        G.build();
        Tree<Graph<T, 1>> tree(G, N);
        return tree;
    }

    // functional graph に向かって進む
    template <typename TREE>
    int jump(TREE& tree, int v, ll step) {
        int d = tree.depth[v];
        if (step <= d - 1) return tree.jump(v, N, step);
        v = root[v];
        step -= d - 1;
        int bottom = TO[v];
        int c = tree.depth[bottom];
        step %= c;
        if (step == 0) return v;
        return tree.jump(bottom, N, step - 1);
    }

    // functional graph に step 回進む
    template <typename TREE>
    vc<int> jump_all(TREE& tree, ll step) {
        auto& G = tree.G;
        vc<int> res(N, -1);
        // v の k 個先を res[w] に入れる
        vvc<pair<int, int>> query(N);
        FOR(v, N) {
            int d = tree.depth[v];
            int r = root[v];
            if (d - 1 > step) { query[v].eb(v, step); }
            if (d - 1 <= step) {
                ll k = step - (d - 1);
                int bottom = TO[r];
                int c = tree.depth[bottom];
                k %= c;
                if (k == 0) {
                    res[v] = r;
                    continue;
                }
                query[bottom].eb(v, k - 1);
            }
        }

        vc<int> path;
        auto dfs = [&](auto & dfs, int v) -> void {
            path.eb(v);
            for (auto && [w, k] : query[v]) { res[w] = path[len(path) - 1 - k]; }
            for (auto && e : G[v]) dfs(dfs, e.to);
            path.pop_back();
        };
        for (auto && e : G[N]) { dfs(dfs, e.to); }
        return res;
    }

    template <typename TREE>
    bool in_cycle(TREE& tree, int v) {
        int r = root[v];
        int bottom = TO[r];
        return tree.in_subtree(bottom, v);
    }

    vc<int> collect_cycle(int r) {
        assert(r == root[r]);
        vc<int> cyc = {TO[r]};
        while (cyc.back() != r) cyc.eb(TO[cyc.back()]);
        return cyc;
    }
};





// Tree Diameter :

vector<int> adj[200001];
int maxx;
int node;
int len;

void dfs(int curr, int parent, int dis) {
    if ( dis > maxx ) {
        node = curr;
        maxx = max(maxx, dis);
    }
    for (auto val : adj[curr]) {
        if ( val != parent)
            dfs(val, curr, dis + 1);
    }
}

void dfs2(int source, int parent, int dis) {
    len = max(len, dis);
    for (auto val : adj[source]) {
        if ( val != parent) {
            dfs2(val, source, dis + 1);
        }
    }
}
void main() {

    int n;
    cin >> n;

    for (int i = 0 ; i < n - 1 ; i++) {
        int a, b;
        cin >> a >> b;
        a--;
        b--;
        adj[a].push_back(b);
        adj[b].push_back(a);
    }

    maxx = -1;
    node = 0;
    dfs(0, -1, 0);
    dfs2(node, -1, 0);

    cout << len << endl; // Length of the tree diameter
}

///