#include "../include/dist_table.hpp"

DistTable::DistTable(const Instance& ins)
    : K(ins.G.V.size()), table(ins.N, std::vector<int>(K, K))
{
  setup(&ins);
}

DistTable::DistTable(const Instance* ins)
    : K(ins->G.V.size()), table(ins->N, std::vector<int>(K, K))
{
  setup(ins);
}

void DistTable::setup(const Instance* ins)
{
  for (size_t i = 0; i < ins->N; ++i) {
    OPEN.push_back(std::queue<Vertex*>());
    auto n = ins->goals[i];
    OPEN[i].push(n);
    table[i][n->id] = 0;
  }
}

int DistTable::get(int i, int v_id)  // 获取agent i到v_id的最短距离
{
  if (table[i][v_id] < K) return table[i][v_id];  // K是顶点数量，表示如果距离被计算过，则直接返回

  /*
   * BFS with lazy evaluation
   * c.f., Reverse Resumable A*
   * https://www.aaai.org/Papers/AIIDE/2005/AIIDE05-020.pdf
   */

  while (!OPEN[i].empty()) {  // i表示agent的编号
    auto n = OPEN[i].front();
    OPEN[i].pop();
    const int d_n = table[i][n->id];
    for (auto& m : n->neighbor) {
      const int d_m = table[i][m->id];
      if (d_n + 1 >= d_m) continue;
      table[i][m->id] = d_n + 1;
      OPEN[i].push(m);
    }
    if (n->id == v_id) return d_n;
  }
  return K;
}

int DistTable::get(int i, Vertex* v) { return get(i, v->id); }
