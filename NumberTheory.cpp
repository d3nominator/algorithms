#include <bits/stdc++.h>
using namespace std;

// prime factorisation with prime and their power
vector<pair<long long, int>> prime_factor(long long x) {
    long long K = x;
    vector<pair<long long, int>> ret;
    for (long long i = 2; i * i <= x; ++i) {
        if (x % i == 0) {
            int cnt = 0;
            while (K % i == 0) {
                K /= i;
                cnt++;
            }
            ret.emplace_back(i, cnt);
        }
    }
    if (K != 1) ret.emplace_back(K, 1);

    return ret;
}


// factors
vector<int> factors(int n) {
    vector<int> f;
    for (int x = 2; x * x <= n; x++) {
        while (n % x == 0) {
            f.push_back(x);
            n /= x;
        }
    }
    if (n > 1) f.push_back(n);
    return f;
}

// Prime factors snippet ^ above is same as this one

vector<int> prime_factor(int n) {
    vector<int> val;
    for (int i = 2 ; i * i <= n ; i++) {
        if ( n % i == 0) {
            while (n % i == 0) {
                val.push_back(i);

                n /= i;
            }
        }
    }
    if (n > 1) {
        val.push_back(n);
    }
    return val;
}

// modular inverse :
inv[1] = 1;
for (int i = 2; i < m; ++i)
    inv[i] = m - (m / i) * inv[m % i] % m;


// segmented sieve :
vector<bool> segmentedSieve(long long L, long long R) {
    // generate all primes up to sqrt(R)
    long long lim = sqrt(R);
    vector<bool> mark(lim + 1, false);
    vector<long long> primes;
    for (long long i = 2; i <= lim; ++i) {
        if (!mark[i]) {
            primes.emplace_back(i);
            for (long long j = i * i; j <= lim; j += i)
                mark[j] = true;
        }
    }

    vector<bool> isPrime(R - L + 1, true);
    for (long long i : primes)
        for (long long j = max(i * i, (L + i - 1) / i * i); j <= R; j += i)
            isPrime[j - L] = false;
    if (L == 1)
        isPrime[0] = false;
    return isPrime;
}

// for computing hash function
long long compute_hash(string const& s) {
    const int p = 31;
    const int m = 1e9 + 9;
    long long hash_value = 0;
    long long p_pow = 1;
    for (char c : s) {
        hash_value = (hash_value + (c - 'a' + 1) * p_pow) % m;
        p_pow = (p_pow * p) % m;
    }
    return hash_value;
}



// binary exponentiation
this is a recursive approaoch
long long binpow(long long a, long long b) {
    if (b == 0)
        return 1;
    long long res = binpow(a, b / 2);
    if (b % 2)
        return res * res * a;
    else
        return res * res;
}

// interative approach for binary exponentiation
ll int binpow(ll int a , ll int b) {
    if (b == 0)
        return 1;
    int res = 1;
    while (b > 0) {
        if (b % 2 == 1)  res = res * a;
        a = a * a;
        b  = b / 2;
    }
    return res;

}

// function to check if a number is a prime number
To find prime numbers
bool isprime(int a) {
    if (a <= 1)
        return false;
    for (int i = 2; i * i < a; i++)
        if (a % i == 0)
            return false;
    return true;
}


//sieve of eratosthenes
** Time complexity : O(n)
for (int x = 2; x <= n; x++) {
    if ( sieve[x] ) continue;
    for (int u = 2 * x; u <= n; u += x) {
        sieve[u] = x;
    }
}



int n;
vector<bool> is_prime(n + 1, true);
is_prime[0] = is_prime[1] = false;
for (int i = 2; i <= n; i++) {
    if (is_prime[i] && (long long)i * i <= n) {
        for (int j = i * i; j <= n; j += i)
            is_prime[j] = false;
    }
}


// number of prime factos a number n had ^-^
int prime_size(int n) {
    int val = 0;
    for (int i = 2 ; i * i <= n ; i++) {
        if ( n % i == 0) {
            while (n % i == 0) {
                val++;
                n /= i;
            }
        }
    }
    if (n > 1) {
        val++;
    }
    return val;
}


// For finding the number of devisors from 1 to n
int n; cin >> n;
vector<int> v(n + 1, 0);
for (int i = 1 ; i <= n ; i++) {
    for (int j = i ; j <= n ; j += i) { // iterate on mulitples instead of devisors
        v[j]++;
    }
}
for (int i = 1 ; i <= n ; i++) {
    cout << v[i] << "\n";
}
// time complexity

/ ------------------------------------------------------------------ /

// Combinatorics :




struct combinatorics
{
    int mod = 1e9 + 7;
    vector<int> fact;
    void init(int n) // necessary
    {
        fact.resize(n + 1);
        fact[0] = 1; fact[1] = 1;
        for (int i = 2; i <= n; i++) fact[i] = ((long long)fact[i - 1] * i) % mod;
    }
    long long powom(long long a, long long b, long long m) {
        a %= m;
        long long res = 1;
        while (b > 0) {
            if (b & 1)
                res = res * a % m;
            a = a * a % m;
            b >>= 1;
        }
        return res;
    }
    int ncrm(int n, int r) // ncr mod m
    {
        int t = (fact[n] * powom(fact[r], mod - 2, mod)) % mod;
        return (t * powom(fact[n - r], mod - 2, mod)) % mod;
    }
    int ncr(int n, int r) { // exact value of ncr
        int ans = 1;
        for (int i = 1; i <= r; i++) {
            ans *= (n + 1 - i);
            ans /= i;
        }
        return ans;
    }
};

/ ---------------------------------------------------------------------------- - /

//Number thoery :
struct number_theory // This template is not mine !!!!
{
    vector<bool>sieve;
    vector<int> prime, spf, etf; // primes , smallest prime factor , euler toient function
    void init(int n)
    {
        sieve.resize(n + 1, 1);
        spf.resize(n + 1);
        etf.resize(n + 1);
        for (int i = 2; i * i <= n; i++) { //O(nlog(log(n)))
            if (sieve[i]) {
                for (int j = i * i; j <= n; j += i) sieve[j] = 0;
            }
        }
        for (int i = 2; i <= n; i++) if (sieve[i]) prime.push_back(i);


        for (int i = 2; i <= n; i++) spf[i] = i;
        for (int i = 2; i * i <= n; i++) { //O(nlog(log(n)))
            if (spf[i] == i) {
                for (int j = i * i; j <= n; j += i) {
                    if (spf[j] == j) spf[j] = i;
                }
            }
        }

        for (int i = 0; i <= n; i++) etf[i] = i;
        for (int i = 2; i <= n; i++) { // O(nlog(n))
            if (etf[i] == i) {
                for (int j = i; j <= n; j += i) {
                    etf[j] -= etf[j] / i;
                }
            }
        }
    }

};


/ ------------------------------------------------------------------------ - /



vector<int> primes;
void pre() {
    vector<bool> is_prime(len + 1, true);
    is_prime[0] = is_prime[1] = false;
    for (int i = 2; i <= len; i++) {
        if (is_prime[i] && (long long)i * i <= len) {
            for (int j = i * i; j <= len; j += i)
                is_prime[j] = false;
        }
    }
    for (int i = 2 ; i <= len ; i ++) {
        if ( is_prime[i]) primes.push_back(i);
    }
}

/ ---------------------------------------------------------------------- - /
//Modular int setup

const int mod = 998244353;
struct mint {
    ll x; // typedef long long ll;
    mint(ll x = 0): x((x % mod + mod) % mod) {}
    mint operator-() const { return mint(-x);}
    mint& operator+=(const mint a) {
        if ((x += a.x) >= mod) x -= mod;
        return *this;
    }
    mint& operator-=(const mint a) {
        if ((x += mod - a.x) >= mod) x -= mod;
        return *this;
    }
    mint& operator*=(const mint a) { (x *= a.x) %= mod; return *this;}
    mint operator+(const mint a) const { return mint(*this) += a;}
    mint operator-(const mint a) const { return mint(*this) -= a;}
    mint operator*(const mint a) const { return mint(*this) *= a;}
    mint pow(ll t) const {
        if (!t) return 1;
        mint a = pow(t >> 1);
        a *= a;
        if (t & 1) a *= *this;
        return a;
    }

    // for prime mod
    mint inv() const { return pow(mod - 2);}
    mint& operator/=(const mint a) { return *this *= a.inv();}
    mint operator/(const mint a) const { return mint(*this) /= a;}
};



struct mi { // WARNING: needs some adjustment to work with FFT
    int v; explicit operator int() const { return v; }
    mi(): v(0) {}
    mi(ll _v): v(int(_v % MOD)) { v += (v < 0) * MOD; }
};
mi& operator+=(mi& a, mi b) {
    if ((a.v += b.v) >= MOD) a.v -= MOD;
    return a;
}
mi& operator-=(mi& a, mi b) {
    if ((a.v -= b.v) < 0) a.v += MOD;
    return a;
}
mi operator+(mi a, mi b) { return a += b; }
mi operator-(mi a, mi b) { return a -= b; }
mi operator*(mi a, mi b) { return mi((ll)a.v * b.v); }
mi& operator*=(mi& a, mi b) { return a = a * b; }
mi pow(mi a, ll p) {
    assert(p >= 0); // won't work for negative p
    return p == 0 ? 1 : pow(a * a, p / 2) * (p & 1 ? a : 1);
}
mi inv(mi a) { assert(a.v != 0); return pow(a, MOD - 2); }
mi operator/(mi a, mi b) { return a * inv(b); }
bool operator==(mi a, mi b) { return a.v == b.v; }




// Custom priority queue with custom delete -- that is for the element not at the top 
template<typename T>
class custom_priority_queue : public std::priority_queue<T, std::vector<T>>
{
  public:

      bool remove(const T& value) {
          auto it = std::find(this->c.begin(), this->c.end(), value);
       
          if (it == this->c.end()) {
              return false;
          }
          if (it == this->c.begin()) {
              // deque the top element
              this->pop();
          }    
          else {
              this->c.erase(it);
              std::make_heap(this->c.begin(), this->c.end(), this->comp);
         }
         return true;
     }
};

