#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <utility>
#include <set>
#include <algorithm>
#include <bits/stdc++.h>

using namespace std;

#define INF 0x3f3f3f3f

class Graph {
    protected:
        int n;
        int m;
        vector<vector<int>> adjList;
        vector<bool> vis;
        vector<int> dist;
    public:
        Graph() {};
};

class UndirectedGraph : public Graph {
    private:
        int noCompConex;
    public:
        UndirectedGraph(string inPath);
        void dfs(int start);
        int dfsCompConex(string outPath);
};

UndirectedGraph :: UndirectedGraph(string inPath) {
    ifstream in(inPath);

    int n, m;

    in >> n >> m;

    this -> n = n;
    this -> m = m;
    adjList.resize(n + 1);
    vis.resize(n + 1);

    for (int i = 0; i < m; i++) {
        int x, y;
        in >> x >> y;
        adjList[x].push_back(y);
        adjList[y].push_back(x);
    }
}

void UndirectedGraph :: dfs(int start) {
    vis[start] = true;

    for (auto x : adjList[start])
        if (!vis[x])
            dfs(x);
}

int UndirectedGraph :: dfsCompConex(string outPath) {
    ofstream out(outPath);

    noCompConex = 0;

    for (int i = 1; i <= n; i++)
        if (!vis[i]) {
            dfs(i);
            noCompConex++;
        }

    out << noCompConex;

    return noCompConex;
}

class DirectedGraph : public Graph {
    private:
        stack<int> topSorted;
    public:
        DirectedGraph(string inPath);
        void bfs(int start);
        vector<int> getBfs(string outPath);
        void dfsTopSort(int start);
        void topSort(string outPath);
};

DirectedGraph :: DirectedGraph(string inPath) {
    ifstream in(inPath);

    int n, m;

    in >> n >> m;

    this -> n = n;
    this -> m = m;
    adjList.resize(n + 1);
    vis.resize(n + 1);
    dist.resize(n + 1, -1);

    for (int i = 0; i < m; i++) {
        int x, y;
        in >> x >> y;
        adjList[x].push_back(y);
    }
}

void DirectedGraph :: bfs(int start) {
    queue<int> q;

    dist[start] = 0;
    q.push(start);

    while (!q.empty()) {
        int node = q.front();
        q.pop();

        for (auto i : adjList[node])
            if (dist[i] == -1) {
                dist[i] = dist[node] + 1;
                q.push(i);
            }
    }
}

vector<int> DirectedGraph :: getBfs(string outPath) {
    ofstream out(outPath);

    for (int i = 1; i <= n; i++)
        out << dist[i] << " ";

    return dist;
}

void DirectedGraph :: dfsTopSort(int start) {
    vis[start] = true;

    for (auto x : adjList[start])
        if (!vis[x])
            dfsTopSort(x);

    topSorted.push(start);
}

void DirectedGraph :: topSort(string outPath) {
    ofstream out(outPath);

    for (int i = 1; i <= n; i++)
        if (!vis[i])
            dfsTopSort(i);

    while (!topSorted.empty()) {
        out << topSorted.top() << " ";
        topSorted.pop();
    }
}

class UndirectedWeightedGraph : public Graph {
    private:
        struct edge {
            int a, b;
            int w;
            friend bool operator<(const edge& A, const edge& B) {
                return A.w < B.w;
            }
        };
        vector<edge> adjList;
        vector<int> subset;
        vector<pair<int, int>> sol;
        int partialCost = 0;
    public:
        UndirectedWeightedGraph(string inPath);
        int findNode(int node);
        void unite(int node1, int node2);
        void kruskal();
        vector<pair<int, int>> apm(string outPath);
        void disjoint(string intPath, string outPath);
};

UndirectedWeightedGraph :: UndirectedWeightedGraph(string inPath) {
    ifstream in(inPath);

    int n, m;

    in >> n >> m;

    this -> n = n;
    this -> m = m;
    subset.resize(n + 1);

    for (int i = 1; i <= n; i++)
        subset[i] = i;

    for (int i = 1; i <= m; i++) {
        edge X;
        in >> X.a >> X.b >> X.w;
        adjList.push_back(X);
    }
}

int UndirectedWeightedGraph :: findNode(int node) {
    if (node == subset[node])
        return node;

    return subset[node] = findNode(subset[node]);
};

void UndirectedWeightedGraph :: unite(int node1, int node2) {
    int a = findNode(node1);
    int b = findNode(node2);
    subset[b] = a;
}

void UndirectedWeightedGraph :: kruskal() {
    sort(adjList.begin(), adjList.end());

    for (auto i : adjList) {
        int a = findNode(i.a);
        int b = findNode(i.b);

        if (a == b)
            continue;
        else {
            unite(i.a, i.b);
            sol.emplace_back(i.a, i.b);
            partialCost += i.w;
        }
    }
}

vector<pair<int, int>> UndirectedWeightedGraph :: apm(string outPath) {
    ofstream out(outPath);

    out << partialCost << '\n' << n - 1 << '\n';
    for (auto i : sol)
        out << i.first << " " << i.second << '\n';

    return sol;
}

void UndirectedWeightedGraph :: disjoint(string inPath, string outPath) {
    ifstream in(inPath);
    ofstream out(outPath);

    int n, m;

    in >> n >> m;

    while (m) {
        int c, x, y;

        in >> c >> x >> y;

        if (c == 1)
            unite(x, y);
        else
            if (findNode(x) == findNode(y))
                out << "DA\n";
            else
                out << "NU\n";
        m--;
    }
}

class DirectedWeightedGraph : public Graph {
    private:
        vector<vector<pair<int, int>>> adjList;
        priority_queue <pair<int, int>> q;
        vector<int> used;
        bool finish = false;
    public:
        DirectedWeightedGraph(string inPath);
        void dijkstra(int start);
        vector<int> dijkstraDist(int start, string outPath);
        void bellmanford(int start);
        vector<int> bellmanfordDist(int start, string outpath);
};

DirectedWeightedGraph :: DirectedWeightedGraph(string inPath) {
    ifstream in(inPath);

    int n, m;

    in >> n >> m;

    this -> n = n;
    this -> m = m;
    adjList.resize(n + 1);
    dist.resize(n + 1, INF);
    used.resize(n + 1);
    vis.resize(n + 1);

    for (int i = 0; i < m; i++) {
        int x, y, w;
        in >> x >> y >> w;
        adjList[x].push_back(make_pair(y, w));
    }
}

void DirectedWeightedGraph :: dijkstra(int start) {
    dist[start] = 0;
    q.push(make_pair(0, start));

    while (!q.empty()) {
        int node = q.top().second;

        if (dist[node] != q.top().first) {
            q.pop();
            continue;
        }

        q.pop();

        for (auto x : adjList[node])
            if (dist[node] + x.second < dist[x.first]) {
                dist[x.first] = dist[node] + x.second;
                q.push(make_pair(dist[x.first], x.first));
            }
    }
}

vector<int> DirectedWeightedGraph :: dijkstraDist(int start, string outPath) {
    ofstream out(outPath);

    for (int i = 1; i <= n; i++)
        if (dist[i] != INF)
            out << dist[i] << " ";
        else
            out << 0 << " ";

    return dist;
}
void DirectedWeightedGraph :: bellmanford(int start) {
    queue<int> q;

    dist[start] = 0;
    q.push(start);

    while (!q.empty()) {
        if (finish)
            break;

        int x = q.front();
        q.pop();
        vis[x] = false;

        for (auto i : adjList[x])
            if (dist[i.first] > dist[x] + i.second) {
                used[i.first]++;

                if (used[i.first]  == n) {
                    finish = true;
                    break;
                }

                dist[i.first] = dist[x] + i.second;

                if (vis[i.first] == false) {
                    q.push(i.first);
                    vis[i.first] = true;
                }
            }
    }
}

vector<int> DirectedWeightedGraph :: bellmanfordDist(int start, string outPath) {
    ofstream out(outPath);

    if (finish)
        out << "Ciclu negativ!";

    for (int i = 1; i <= n; i++)
        out << dist[i] << " ";

    return dist;
}

class Tree : public Graph {
    private:
        int diam;
        int last;
    public:
        Tree(string inPath);
        void bfsDarb(int start);
        int distDarb(int start, string outPath);
};

Tree :: Tree(string inPath) {
    ifstream in(inPath);

    int n;

    in >> n;

    this -> n = n;
    adjList.resize(n + 1);
    dist.resize(n + 1);

    for (int i = 1; i < n; i++) {
        int x, y;
        in >> x >> y;
        adjList[x].push_back(y);
        adjList[y].push_back(x);
    }
}

void Tree :: bfsDarb(int start) {
    queue<int> q;
    dist[start] = 1;
    q.push(start);

    while (!q.empty()) {
        int x = q.front();
        last = x;

        for (auto y : adjList[x]) {
            if (!dist[y]) {
                dist[y] = dist[x] + 1;
                q.push(y);
            }
        }

        q.pop();
    }

    diam = dist[last];

}

int Tree :: distDarb(int start, string outPath) {
    ofstream out(outPath);

    bfsDarb(start);
    dist.assign(n + 1, 0);
    bfsDarb(last);

    out << diam;

    return diam;
}

class RoyFloydGraph : public Graph {
    public:
        RoyFloydGraph(string inPath);
        void royFloyd();
        vector<vector<int>> getRoyFloyd(string outPath);

};

RoyFloydGraph :: RoyFloydGraph(string inPath) {
    ifstream in(inPath);

    int n;

    in >> n;

    this -> n = n;
    adjList.resize(n + 1, vector<int> (n + 1, 0));

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            in >> adjList[i][j];
            if (adjList[i][j] == 0 && i != j)
                adjList[i][j] = INF;
        }
}

void RoyFloydGraph :: royFloyd() {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            for (int k = 0; k < n; k++)
                if (adjList[j][i] + adjList[i][k] < adjList[j][k])
                    adjList[j][k] = adjList[j][i] + adjList[i][k];
}

vector<vector<int>> RoyFloydGraph :: getRoyFloyd(string outPath) {
    ofstream out(outPath);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            out << adjList[i][j] << " ";
        out << '\n';
    }

    return adjList;
}

class ResidualGraph : public Graph {
    private:
        vector<vector<int>> capacity;
        vector<vector<int>> flow;
        vector<int> pred;
        int flowValue = 0;
    public:
        ResidualGraph(string inPath);
        int getN();
        bool bfsResidualGraph(int start, int dest);
        int findMaxFlow(int start, int dest, string outPath);
};

ResidualGraph :: ResidualGraph(string inPath) {
    ifstream in(inPath);

    int n, m;

    in >> n >> m;

    this -> n = n;
    this -> m = m;
    adjList.resize(n + 1);
    vis.resize(n + 1);
    pred.resize(n + 1);
    capacity.resize(n + 1, vector<int> (n + 1, 0));
    flow.resize(n + 1, vector<int> (n + 1, 0));

    for (int i = 0; i < m; i++) {
        int x, y, c;
        in >> x >> y >> c;
        adjList[x].push_back(y);
        adjList[y].push_back(x);
        capacity[x][y] = c;
    }
}

int ResidualGraph :: getN() {
    return n;
}

bool ResidualGraph :: bfsResidualGraph(int start, int dest) {
    queue<int> q;

    vis.assign(n + 1, 0);

    q.push(start);
    vis[start] = true;

    while (!q.empty()) {
        int node = q.front();
        q.pop();

        if (node == dest)
            continue;

        for (auto x : adjList[node]) {
            if (capacity[node][x] == flow[node][x])
                continue;

            if (vis[x])
                continue;

            q.push(x);
            pred[x] = node;
            vis[x] = true;
        }
    }

    return vis[dest];
}

int ResidualGraph :: findMaxFlow(int start, int dest, string outPath) {
    ofstream out(outPath);

    while (bfsResidualGraph(start, dest))
        for (auto x : adjList[dest]) {
            int flowMax = INF;

            if (capacity[x][dest] == flow[x][dest])
                continue;

            if (!vis[x])
                continue;

            pred[dest] = x;

            for (int node = dest; node != start; node = pred[node])
                flowMax = min(flowMax, capacity[pred[node]][node] - flow[pred[node]][node]);

            if (!flowMax)
                continue;

            for (int node = dest; node != start; node = pred[node]) {
                flow[pred[node]][node] += flowMax;
                flow[node][pred[node]] -= flowMax;
            }

            flowValue += flowMax;
        }

    out << flowValue;

    return flowValue;
}
///////////////////////////////

int main()
{
    /* dfs O(n + m)
    UndirectedGraph X = UndirectedGraph("dfs.in");
    X.dfsCompConex("dfs.out");
    */

    /* bfs O(n + m)
    DirectedGraph X = DirectedGraph("bfs.in");
    X.bfs(2);
    X.getBfs("bfs.out");
    */

    /* Top Sort O(n + m)
    DirectedGraph X = DirectedGraph("sortaret.in");
    X.topSort("sortaret.out");
    */

    /// Tema 2

    /* apm O(m*logn + m*logm)
    UndirectedWeightedGraph X = UndirectedWeightedGraph("apm.in");
    X.kruskal();
    X.apm("apm.out");
    */

    /* disjoint O(m)
    UndirectedWeightedGraph X("disjoint.in");
    X.disjoint("disjoint.in", "disjoint.out");
    */

    /* dijkstra O(m*logn)
    DirectedWeightedGraph X = DirectedWeightedGraph("dijkstra.in");
    X.dijkstra(1);
    X.dijkstraDist(1, "dijkstra.out");
    */

    /* bellman ford O(n*m)
    DirectedWeightedGraph X = DirectedWeightedGraph("bellmanford.in");
    X.bellmanford(1);
    X.bellmanfordDist(1, "bellmanford.out");
    */

    /// Tema 3

    /*Flow Max O(n*m^2)
    ResidualGraph X = ResidualGraph("maxflow.in");
    X.findMaxFlow(1, X.getN(), "maxflow.out");
    */

    /* darb O(n + m)
    Tree X = Tree("darb.in");
    X.distDarb(1, "darb.out");
    */
    /* Roy Floyd O(n^3)
    RoyFloydGraph X = RoyFloydGraph("royfloyd.in");
    X.royFloyd();
    X.getRoyFloyd("royfloyd.out");
    */

    return 0;
}
