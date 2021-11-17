#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <utility>
#include <set>
#include <algorithm>

using namespace std;

#define INF 0x3f3f3f3f

class Graph {
    private:
        int n;
        int m;
        int start;
        int c;
        static int noBC;
        vector<vector<int>> adjList;
        vector<vector<int>> adjListT;
        vector<int> dist;
        queue<int> q;
        bool vis[100001] = {false};
        bool visT[100001] = {false};
        vector<int> solComp[100001];
        stack<int> S;
        vector<int> lvl;
        vector<int> low;
        vector<vector<int>> solBC;
        vector<int> subset;
    public:
        Graph() {};
        void directedGraph();
        void undirectedGraph();
        void undirectedGraphBC();
        void directedGraphTopSort();
        void directedGraphCTC();
        void bfs();
        void dfs(int start);
        void noConnComp();
        void dfsBC(int start, int pred);
        int getSolBC();
        void dfsTopSort(int start, stack<int>& topSorted);
        void topSort();
        void dfsT(int start);
        void printCTC();
        void dijkstra();
        void bellmanford();
        struct edge {
            int a, b;
            int w;
            friend bool operator<(const edge& A, const edge& B) {
                return A.w < B.w;
            }
        };
        int findNode(int node);
        void unite(int node1, int node2);
        void apm();
        void disjoint();
};

int Graph :: noBC = 0;

void Graph :: directedGraph() {
    ifstream fin("bfs.in");

    fin >> n >> m >> start;

    adjList.resize(n + 1);

    for (int i = 0; i < m; i++) {
        int x, y;
        fin >> x >> y;
        adjList[x].push_back(y);
    }

    fin.close();
}

void Graph :: undirectedGraph() {
    ifstream fin("dfs.in");

    fin >> n >> m;

    adjList.resize(n + 1);

    for (int i = 0; i < m; i++) {
        int x, y;
        fin >> x >> y;
        adjList[x].push_back(y);
        adjList[y].push_back(x);
    }

    fin.close();
}

void Graph :: undirectedGraphBC() {
    ifstream fin("biconex.in");

    fin >> n >> m;

    adjList.resize(n + 1);
    solBC.resize(n + 1);
    lvl.resize(n + 1);
    low.resize(n + 1);

    for (int i = 0; i < m; i++) {
        int x, y;
        fin >> x >> y;
        adjList[x].push_back(y);
        adjList[y].push_back(x);
    }

    fin.close();
}

void Graph :: directedGraphTopSort() {
    ifstream fin("sortaret.in");

    fin >> n >> m;

    adjList.resize(n + 1);

    for (int i = 0; i < m; i++) {
        int x, y;
        fin >> x >> y;
        adjList[x].push_back(y);
    }

    fin.close();
}

void Graph :: directedGraphCTC() {
    ifstream fin("ctc.in");

    fin >> n >> m;

    adjList.resize(n + 1);
    adjListT.resize(n + 1);

    for (int i = 0; i < m; i++) {
        int x, y;
        fin >> x >> y;
        adjList[x].push_back(y);
        adjListT[y].push_back(x);
    }

    fin.close();
}

void Graph :: bfs() {
    ofstream fout("bfs.out");

    q.push(start);

    dist.resize(n + 1);
    for (int i = 1; i <= n; i++)
        dist[i] = -1;
    dist[start] = 0;

    while (!q.empty()) {
        int nod = q.front();
        q.pop();
        for (auto i : adjList[nod])
            if (dist[i] == -1) {
                dist[i] = dist[nod] + 1;
                q.push(i);
            }
    }

    for (int i = 1; i <= n; i++)
        fout << dist[i] << " ";

    fout.close();
}

void Graph :: dfs(int start) {
    vis[start] = true;

    for (auto x : adjList[start])
        if (!vis[x])
            dfs(x);
}

void Graph :: noConnComp() {
    ofstream fout("dfs.out");

    int c = 0;

    for (int i = 1; i <= n; i++)
        if (!vis[i]) {
            dfs(i);
            c++;
        }

    fout << c;

    fout.close();
}

void Graph :: dfsBC(int start, int pred) {
    vis[start] = 1;
    S.push(start);
    lvl[start] = low[start] = lvl[pred] + 1;

    for (auto x : adjList[start])
        if (x != pred)
            if (vis[x])
                if (low[start] > lvl[x])
                    low[start] = lvl[x];
            else {
                dfsBC(x, start);

                if (low[start] > low[x])
                    low[start] = low[x];

                if (lvl[start] < low[x]) {
                    noBC++;

                    while (S.top() != x) {
                        solBC[noBC].push_back(S.top());
                        S.pop();
                    }

                    solBC[noBC].push_back(x);
                    S.pop();
                    solBC[noBC].push_back(start);
                }
            }

}

int Graph :: getSolBC() {
    return noBC;
}

void Graph :: dfsTopSort(int start, stack<int>& topSorted) {
    vis[start] = true;

    for (auto x : adjList[start])
        if (!vis[x])
            dfsTopSort(x, topSorted);

    topSorted.push(start);
}

void Graph :: topSort() {
    ofstream fout("sortaret.out");

    stack<int> topSorted;

    for (int i = 1; i <= n; i++)
        if (!vis[i])
            dfsTopSort(i, topSorted);

    while (!topSorted.empty()){
        fout << topSorted.top() << " ";
        topSorted.pop();
    }

    fout.close();
}

void Graph :: dfsT(int start) {
    visT[start] = true;

    solComp[c].push_back(start);

    for (auto x : adjList[start])
        if (!visT[x])
            dfsT(x);
}

void Graph :: printCTC() {
    ofstream fout("ctc.out");

    stack<int> topSorted;

    for (int i = 1; i <= n; i++)
        if (!vis[i])
            dfsTopSort(i, topSorted);

    while (!topSorted.empty()) {
        if (!visT[topSorted.top()]) {
            c++;
            dfsT(topSorted.top());
        }
        topSorted.pop();
    }

    cout << c << '\n';

    for (int i = 0; i < c; i++) {
        for (int j = 0; j < solComp[i].size(); j++)
            cout << solComp[i][j] << " ";
        cout << '\n';
    }

    fout.close();
}

void Graph :: dijkstra() {
    ifstream fin("dijkstra.in");
    ofstream fout("dijkstra.out");

    vector<list<pair<int, int>>> adjList;
    int n, m;

    fin >> n >> m;

    adjList.resize(n + 1);

    for (int i = 0; i < m; i++) {
        int a, b, w;
        fin >> a >> b >> w;
        adjList[a].push_back(make_pair(b, w));
    }

    set<pair<int, int>> setProc;
    vector<int> dist(n + 1, INF);

    dist[1] = 0;
    setProc.insert(make_pair(0, 1));

    while (!setProc.empty()) {
        pair<int, int> tmp = *(setProc.begin());
        setProc.erase(setProc.begin());

        int u = tmp.second;

        list<pair<int, int>> :: iterator i;
        for (i = adjList[u].begin(); i != adjList[u].end(); i++) {
            int v = (*i).first;
            int w = (*i).second;

            if (dist[v] > dist[u] + w) {
                if (dist[v] != INF)
                    setProc.erase(setProc.find(make_pair(dist[v], v)));

                dist[v] = dist[u] + w;
                setProc.insert(make_pair(dist[v], v));
            }
        }
    }

    for (int i = 2; i <= n; i++)
        fout << dist[i] << " ";
}

void Graph :: bellmanford() {
    ifstream fin("bellmanford.in");
    ofstream fout("bellmanford.out");

    int n, m;
    vector<vector<pair<int, int>>> adjList;

    fin >> n >> m;

    adjList.resize(n + 1);

    for (int i = 0; i < m; i++) {
        int a, b, w;
        fin >> a >> b >> w;
        adjList[a].push_back(make_pair(b, w));
    }

    vector<int> dist(n + 1, INF);
    queue<int> q;
    vector<bool> vis;
    vector<int> used;

    vis.resize(n + 1);

    used.resize(n + 1);

    dist[1] = 0;
    q.push(1);

    bool finish = false;

    while (!q.empty()) {
        if (finish)
            break;

        int x = q.front();
        q.pop();
        vis[x] = false;

        for (auto i : adjList[x])
            if (dist[i.first] > dist[x] + i.second) {
                ++used[i.first];

                if (used[i.first]  == n) {
                    fout << "Ciclu negativ!";
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

    if (!finish)
        for (int i = 2; i <= n; i++)
            fout << dist[i] << " ";
}

int Graph :: findNode(int node) {
    if (node == subset[node])
        return node;
    return subset[node] = findNode(subset[node]);
};

void Graph :: unite(int node1, int node2) {
    int a = findNode(node1);
    int b = findNode(node2);
    subset[b] = a;
}

void Graph :: apm() {
    ifstream fin("apm.in");
    ofstream fout("apm.out");

    int n, m;
    vector<edge> adjList;
    vector<pair<int, int>> sol;

    fin >> n >> m;

    subset.resize(n + 1);

    for (int i = 1; i <= n; i++)
        subset[i] = i;

    for (int i = 1; i <= m; i++) {
        edge X;
        fin >> X.a >> X.b >> X.w;
        adjList.push_back(X);
    }

    sort(adjList.begin(), adjList.end());

    int s = 0;

    for (auto i : adjList) {
        int a = findNode(i.a);
        int b = findNode(i.b);

        if (a == b)
            continue;
        else {
            unite(i.a, i.b);
            sol.emplace_back(i.a, i.b);
            s += i.w;
        }
    }

    fout << s << '\n' << n - 1 << '\n';
    for (auto i : sol)
        fout << i.first << " " << i.second << '\n';

}

void Graph :: disjoint() {
    ifstream fin("disjoint.in");
    ofstream fout("disjoint.out");

    int n, m;

    fin >> n >> m;

    subset.resize(n + 1);

    for (int i = 1; i <= n; i++)
        subset[i] = i;

    while (m) {
        int c, x, y;

        fin >> c >> x >> y;

        if (c == 1)
            unite(x, y);
        else
            if (findNode(x) == findNode(y))
                fout << "DA\n";
            else
                fout << "NU\n";
        m--;
    }
}

int main()
{
    Graph X;
    X.disjoint();
    return 0;
}
