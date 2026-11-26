// Enumerate small, irreducible detector-silent supports with nonzero payload.
// Input: detector width, sorted column count, then unsigned low/high limbs.
// Output: each support as column indices, followed by the number of visited nodes.
#include <algorithm>
#include <bitset>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

using Word = unsigned __int128;
constexpr int capacity = 512;
using Sites = std::bitset<capacity>;

struct Hash {
    std::size_t operator()(Word value) const {
        return std::uint64_t(value) ^ (std::uint64_t(value >> 64) * 11400714819323198485ull);
    }
};

class Search {
public:
    Search(int width, const std::vector<Word>& columns, std::uint64_t budget)
        : cols(columns), det((Word(1) << width) - 1), checks(width), sites(width), budget(budget) {
        for (int i = 0; i < int(cols.size()); ++i) {
            for (int c = 0; c < width; ++c) {
                if ((cols[i] >> c) & 1) {
                    checks[c].set(i);
                    sites[c].push_back(i);
                }
            }
            tails[cols[i] & det].push_back({i, -1});
            for (int j = 0; j < i; ++j)
                tails[(cols[i] ^ cols[j]) & det].push_back({j, i});
        }
    }

    void run(int limit) {
        Sites free;
        for (int q = 0; q < int(cols.size()); ++q) free.set(q);
        for (int q = 0; q < int(cols.size()); ++q) {
            free.reset(q);
            prefix = {q};
            visit(cols[q], free, limit - 1);
        }
        std::cout << "nodes " << visited << '\n';
    }

private:
    const std::vector<Word>& cols;
    Word det;
    std::vector<Sites> checks;
    std::vector<std::vector<int>> sites;
    std::unordered_map<Word, std::vector<std::pair<int,int>>, Hash> tails;
    std::uint64_t budget, visited = 0;
    std::vector<int> prefix;

    void emit(int a = -1, int b = -1) const {
        for (int i : prefix) std::cout << i << ' ';
        if (a >= 0) std::cout << a << ' ';
        if (b >= 0) std::cout << b;
        std::cout << '\n';
    }

    void visit(Word sum, Sites free, int left) {
        if (++visited > budget) throw std::runtime_error("node budget");
        Word odd = sum & det;
        if (!odd) {
            if (sum) emit();
            return;
        }
        if (left == 2) {
            auto found = tails.find(odd);
            if (found == tails.end()) return;
            for (auto [a,b] : found->second) {
                if (!free[a] || (b >= 0 && !free[b])) continue;
                if (sum ^ cols[a] ^ (b >= 0 ? cols[b] : 0)) emit(a,b);
            }
            return;
        }
        int pivot = -1, smallest = capacity+1, packed = 0;
        Sites used;
        for (int offset : {0,64}) {
            std::uint64_t limb = std::uint64_t(odd >> offset);
            while (limb) {
                int c = offset + __builtin_ctzll(limb);
                limb &= limb - 1;
                Sites options = free & checks[c];
                int count = int(options.count());
                if (count == 0) return;
                if (count < smallest) { smallest = count; pivot = c; }
                if ((options & used).none()) {
                    if (++packed > left) return;
                    used |= options;
                }
            }
        }
        // Partition by the first selected available column in the chosen row.
        for (int i : sites[pivot]) {
            if (!free[i]) continue;
            free.reset(i);
            prefix.push_back(i);
            visit(sum ^ cols[i], free, left-1);
            prefix.pop_back();
        }
    }
};

int main(int argc, char** argv) {
    int width, n;
    if (argc != 3 || !(std::cin >> width >> n) || width < 1 || width >= 128 ||
        n < 1 || n > capacity) return 2;
    int limit = std::stoi(argv[1]);
    if (limit < 3 || limit > 32) return 2;
    std::vector<Word> cols;
    for (int i = 0; i < n; ++i) {
        std::uint64_t lo,hi;
        if (!(std::cin >> lo >> hi)) return 2;
        cols.push_back((Word(hi) << 64) | lo);
    }
    try {
        Search(width,cols,std::stoull(argv[2])).run(limit);
    } catch (const std::runtime_error&) {
        return 3;
    }
}
