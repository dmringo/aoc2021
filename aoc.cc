#include <algorithm>
#include <bitset>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <ostream>
#include <regex>
#include <set>
#include <sstream>
#include <tuple>
#include <variant>
#include <vector>

#define LAST_DAY 25
#if DEBUG
#define DBG(x) x
#else
#define DBG(x)
#endif

using std::string;

using AOCResult = std::variant<long, string>;
using AOCSoln = std::function<AOCResult(std::istream&)>;

/// Run a function for every line in the input istream. Function can accumulate
/// a value in its second parameter, and should return false if iteration should
/// stop, true otherwise.
template<typename Acc_t>
void
maplines(std::istream& file,
         std::function<bool(string&, Acc_t&)> func,
         Acc_t& acc)
{
  string line;
  while (std::getline(file, line) && func(line, acc)) {
    // pass
  }
}

/// Wrapper of maplines above without accumulator
void
maplines(std::istream& file, std::function<bool(string&)> func)
{
  bool _ignored = false;
  auto wrapper_func = [&](string l, bool&) { return func(l); };
  maplines<bool>(file, wrapper_func, _ignored);
}

/// ostream operator<< for variants of types that have an ostream operator<<
/// impl
// Shamelessly stolen from https://stackoverflow.com/a/47168851/8541141
template<typename T0, typename... Ts>
std::ostream&
operator<<(std::ostream& s, std::variant<T0, Ts...> const& v)
{
  std::visit([&](auto&& arg) { s << arg; }, v);
  return s;
}

/// Day 1 Part 1
///
/// Input file has a single number of each line, representing an ocean depth
/// reading. Goal is to simply count the number of measurements that increase in
/// depth from the prior. For the first line, there is no prior measurement.
AOCResult
day1p1(std::istream& input)
{
  struct Acc
  {
    /// Count of increasing measurements
    unsigned total = 0;
    // whatever the last measurement was
    unsigned last = UINT_MAX;
  } acc;

  auto count_increasing = [&acc](string& line) {
    unsigned current = std::stoi(line);
    DBG(std::cerr << "Last: " << acc.last << ", Current: " << current << "\n");
    acc.total += current > acc.last;
    acc.last = current;
    return true;
  };
  maplines(input, count_increasing);

  return acc.total;
};

/// Day 1 Part 2
///
/// Input format same as Day 1 Part 1. Now the challenge is to find when the
/// rolling sum of the last three lines increases.
AOCResult
day1p2(std::istream& input)
{
  /// Count of increasing sums
  unsigned total = 0;
  /// Prior depths
  unsigned depths[3];
  /// Line number
  unsigned linum = 0;

  auto count_increasing = [&](string& line) {
    if (linum > 2) {
      unsigned last = depths[0] + depths[1] + depths[2];
      depths[linum % 3] = std::stoi(line);
      unsigned current = depths[0] + depths[1] + depths[2];
      total += current > last;
    } else {
      depths[linum % 3] = std::stoi(line);
    }
    linum++;
    return true;
  };

  maplines(input, count_increasing);

  return total;
}

enum class Dir
{
  Forward,
  Up,
  Down,
  Bad
};

static auto
map_day_2(std::istream& input, std::function<bool(Dir, unsigned)> func)
{
  static const std::regex rx{ R"((\S+) (\d+))",
                              std::regex_constants::ECMAScript };
  std::smatch match;
  maplines(input, [&](string line) {
    if (std::regex_match(line, match, rx)) {
      string dir = match[1];
      Dir edir = Dir::Bad;
      unsigned mag = std::stoi(match[2]);
      if (dir == "forward") {
        edir = Dir::Forward;
      } else if (dir == "up") {
        edir = Dir::Up;
      } else if (dir == "down") {
        edir = Dir::Down;
      } else {
        DBG(std::cerr << "Unrecognized direction: " << std::quoted(dir)
                      << "\n");
      }
      return func(edir, mag);
    } else {
      DBG(std::cerr << "Unmatched line???\n" << line << "\n");
      return false;
    }
  });
}
/** Day 2 Part 1
 *
 * Input is lines with a direction {forward, down, up} and a number. 'forward'
 * increases horizontal position, down increases depth, up decreases depth.
 * Calculate final position (assuming (0,0) starting position) and return the
 * *product* of the two components.
 */
AOCResult
day2p1(std::istream& input)
{
  unsigned depth = 0;
  unsigned hpos = 0;
  map_day_2(input, [&](Dir d, unsigned mag) {
    switch (d) {
      case Dir::Forward:
        hpos += mag;
        break;
      case Dir::Up:
        depth -= mag;
        break;
      case Dir::Down:
        depth += mag;
        break;
      case Dir::Bad:
        DBG(std::cerr << "BAD??\n");
        return false;
    }
    return true;
  });

  return depth * hpos;
}

/** Day 2 Part 2
 *
 * Now an 'up' or 'down' directive influences the submarines aim, instead of
 * depth directly. Aim is increased by X with 'down' and decreased by X with
 * 'up'. When 'forward' is encountered, horizontal position is increased by X,
 * *and* depth is increased by X * aim. We still return the product of
 * horizontal position and depth as the answer.
 * */
AOCResult
day2p2(std::istream& input)
{
  unsigned aim = 0;
  unsigned depth = 0;
  unsigned hpos = 0;
  map_day_2(input, [&](Dir d, unsigned mag) {
    switch (d) {
      case Dir::Forward:
        hpos += mag;
        depth += mag * aim;
        break;
      case Dir::Up:
        aim -= mag;
        break;
      case Dir::Down:
        aim += mag;
        break;
      case Dir::Bad:
        DBG(std::cerr << "BAD??\n");
        return false;
    }
    return true;
  });
  return depth * hpos;
}

#define NBITS 12
/** Day 3 Part 1
 *
 * Input lines are binary numbers (on inspection, 12 bits). Goal is to produce
 * two numbers (gamma and epsilon) that are composed of the most and least
 * common bits at each position across all lines, respectively. Return the
 * product of these two numbers as the answer.
 */
AOCResult
day3p1(std::istream& input)
{
  // Assumptions:
  // - Lines all have exactly 12 characters, all of which are either '1' or '0'
  // - Based on question phrasing, there won't be bits at which 1 and 0 are
  //   equally frequent.
  // - There aren't more than UINT_MAX lines in the input

  /// Store count of 1s seen for each bit index. Left most bit is highest index
  unsigned bit_counts[NBITS] = { 0 };
  /// count number of lines to compare for most frequent
  unsigned nlines = 0;

  maplines(input, [&](string& line) {
    nlines++;
    unsigned i = NBITS - 1;
    for (auto& ch : line) {
      bit_counts[i--] += ch == '1';
    }
    return true;
  });

  unsigned gamma = 0;
  unsigned epsilon = 0;
  for (int i = NBITS - 1; i >= 0; i--) {
    // 1 is most common i-th bit, set in gamma
    DBG(std::cerr << "Bit " << std::setw(2) << i << " : " << bit_counts[i]
                  << "\n");

    if (bit_counts[i] > (nlines / 2)) {
      gamma |= 1 << i;
    } else {
      epsilon |= 1 << i;
    }
  }
#ifdef DEBUG
  std::bitset<12> g(gamma);
  std::bitset<12> e(epsilon);
  std::cerr << "Gamma: " << g << "\n"
            << "Epsilon: " << e << "\n";
#endif

  return gamma * epsilon;
}

/** Day 3 Part 2
 *
 * Now two numbers (oxygen generator rating, CO2 scrubber rating) are to be
 * found, based on "bit criteria". Starting with the entire set of numbers, keep
 * only those with the most common value in the leftmost bit in the set. If both
 * values occur with the same frequency, keep those with a 1. If one number
 * remains, that's the oxygen generator rating. Otherwise, repeat with the next
 * bit (considering only frequencies in the remainder of the set), and so on.
 *
 * CO2 scrubber rating is the same process, but keeping numbers with the *least*
 * common value, and breaking a tie with numbers holding 0s.
 *
 * The answer is the product of the two resulting numbers.
 */
AOCResult
day3p2(std::istream& input)
{
  // We probably need O(nlines) space, but can do better than O(n) time for the
  // successive filtering steps by taking advantage of binary search.
  // Tree structure represents bits in a number. Insertion of a number follows
  // bits from most to least significant.
  struct Tree
  {
    std::unique_ptr<Tree> children[2];
    unsigned count = 0;

    // Insert number n into tree. idx is bit index: should be NBITS-1 when
    // called at root
    void insert(unsigned n, int idx = NBITS - 1)
    {
      count++;
      if (idx < 0)
        return;
      unsigned bit = (n >> idx) & 1;
      if (nullptr == children[bit]) {
        children[bit] = std::make_unique<Tree>();
      }
      children[bit]->insert(n, idx - 1);
    }

    // Find number given a comparator on the children counts
    unsigned find(std::function<bool(unsigned, unsigned)> cmp,
                  int idx = NBITS - 1)
    {
      unsigned bit = 0;
      if (children[0] == nullptr) {
        if (children[1] == nullptr) {
          // base case, both are absent
          // should assert idx == -1, probably
          return 0;
        }
        // child 0 null => only child 1 exists
        bit = 1;
      } else if (children[1] != nullptr) {
        // both are present, have to compare counts
        bit = cmp(children[0]->count, children[1]->count);
      }

      DBG(std::cerr << "find(idx = " << idx << ") - count = " << std::setw(8)
                    << count << "\n  "
                    << "bit = " << bit << "\n");

      return children[bit]->find(cmp, idx - 1) | (bit << idx);
    }
  };

  auto root = std::make_unique<Tree>();

  maplines(input, [&](string& line) {
    unsigned num = 0;
    int i = NBITS - 1;
    for (auto ch : line) {
      unsigned bit = ch == '1';
      num |= bit << i;
      i--;
    }
    root->insert(num);
    return true;
  });

  unsigned ox_gen_rating =
    root->find([&](unsigned c0, unsigned c1) { return c1 >= c0; });
  unsigned co2_scrub_rating =
    root->find([&](unsigned c0, unsigned c1) { return c1 < c0; });

  return ox_gen_rating * co2_scrub_rating;
}

#undef NBITS

#define NROWS 5
#define COL_START NROWS
#define DIAG_START (COL_START + NROWS)

struct Board
{
  using Line = std::set<unsigned>;
  // lines in board represented as sets of unsigneds.
  // First 5 are the rows.
  // Following 5 are the columns
  // Last two are the negative and positive sloped diagonals, in order.
  Line lines[12] = { {} };

  // keep track of whether this board has won - important for part 2
  bool bingo = false;

  Line& operator[](int i) { return lines[i]; }

  void add_row(string s, int row)
  {
    size_t pos;
    for (int col = 0; col < NROWS; col++) {
      // parse number, get position
      unsigned n = std::stoi(s, &pos);
      // erase digits and whatever space follows
      s.erase(0, pos + 1);
      lines[row].insert(n);
      lines[COL_START + col].insert(n);
      if (col == row) {
        lines[DIAG_START].insert(n);
      }
      if (row + col == NROWS - 1) {
        lines[DIAG_START + 1].insert(n);
      }
    }
  }

  void fill(std::istream& input)
  {
    string s;
    for (int row = 0; row < NROWS; row++) {
      std::getline(input, s);
      add_row(s, row);
    }
  }

  // Helper for reading first line of input, doesn't belong here, but that's OK
  static std::vector<unsigned> read_called_nums(std::istream& input)
  {

    std::vector<unsigned> nums;
    string s;
    size_t pos;
    std::getline(input, s);
    do {

      // read number (comma/end of string will halt parsing in stoi)
      nums.push_back(std::stoi(s, &pos));
      // remove number *and* comma
      s.erase(0, pos + 1);
      // s.empty() when last number removed from s in final iteration
    } while (!s.empty());

    return nums;
  }
};
/** Day 4 Part 1
 *
 * Bingo! Input is a line of comma-delimited numbers, followed by several
 * occurences of a blank line followed by a 5x5 space-and-newline-delimited
 * board.
 *
 * The numbers in the first line are the numbers called, in order, for a Bingo
 * game. Find the first board that will get a Bingo. The answer is the sum of
 * the *unmarked* spots on that board multiplied by the final number called.
 */
AOCResult
day4p1(std::istream& input)
{
  // first, parse the line of called numbers
  std::vector<unsigned> nums = Board::read_called_nums(input);

  // Discard newline - should probably assert ret == '\n'
  input.get();

  std::vector<Board> boards;
  while (!input.eof()) {
    Board b;
    b.fill(input);

    boards.push_back(b);
    // Discard newline - should probably assert ret == '\n'
    input.get();
  }

  // Now play Bingo
  Board winner;
  unsigned winning_num;
  bool done = false;
  for (auto n : nums) {
    for (auto& board : boards) {
      for (auto& line : board.lines) {
        line.erase(n);
        if (line.size() == 0) {
          done = true;
          winner = board;
          winning_num = n;
        }
        // keep going to erase the number from the remaining lines
      }
      if (done)
        break;
    }
    if (done)
      break;
  }

  std::set<unsigned> unmarked;
  for (auto& line : winner.lines) {
    unmarked.merge(line);
  }
  unsigned sum = 0;
  for (auto n : unmarked) {
    sum += n;
  }

  return sum * winning_num;
}

/** Day 4 Part 2
 *
 * Now we want to find the *last* board to win. As before, the answer is the sum
 * of the unmarked spots on that board multiplied by the last number
 */
AOCResult
day4p2(std::istream& input)
{
  // first, parse the line of called numbers
  std::vector<unsigned> nums = Board::read_called_nums(input);

  // Discard newline - should probably assert ret == '\n'
  input.get();

  std::vector<Board> boards;
  while (!input.eof()) {
    Board b;
    b.fill(input);

    boards.push_back(b);
    // Discard newline - should probably assert ret == '\n'
    input.get();
  }

  unsigned n, ni = 0;
  while (boards.size() > 1) {
    n = nums[ni++];
    for (auto& board : boards) {
      for (auto& line : board.lines) {
        line.erase(n);
        if (line.size() == 0) {
          board.bingo = true;
          break;
        }
      }
    }
    boards.erase(std::remove_if(boards.begin(),
                                boards.end(),
                                [](auto board) { return board.bingo; }),
                 boards.end());
  }

  auto& winner = boards[0];
  while (!winner.bingo) {
    n = nums[ni++];
    for (auto& line : winner.lines) {
      line.erase(n);
      if (line.size() == 0) {
        winner.bingo = true;
        // keep going, to erase n from any other lines it might be in
      }
    }
  }

  std::set<unsigned> unmarked;
  unsigned sum = 0;
  const unsigned winning_num = n;
  for (auto& line : winner.lines) {
    for(auto spot : line) {
      if (unmarked.insert(spot).second) {
        sum += spot;
      }
    }
  }
  return sum * winning_num;
}

AOCResult
day5p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day5p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day6p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day6p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day7p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day7p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day8p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day8p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day9p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day9p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day10p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day10p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day11p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day11p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day12p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day12p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day13p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day13p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day14p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day14p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day15p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day15p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day16p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day16p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day17p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day17p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day18p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day18p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day19p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day19p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day20p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day20p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day21p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day21p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day22p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day22p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day23p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day23p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day24p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day24p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day25p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day25p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day26p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day26p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day27p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day27p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day28p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day28p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day29p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day29p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day30p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day30p2(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day31p1(std::istream& input)
{
  return "Unimplemented!";
}

AOCResult
day31p2(std::istream& input)
{
  return "Unimplemented!";
}

static std::vector<std::vector<AOCSoln>> days = {
  { day1p1, day1p2 },   { day2p1, day2p2 },   { day3p1, day3p2 },
  { day4p1, day4p2 },   { day5p1, day5p2 },   { day6p1, day6p2 },
  { day7p1, day7p2 },   { day8p1, day8p2 },   { day9p1, day9p2 },
  { day10p1, day10p2 }, { day11p1, day11p2 }, { day12p1, day12p2 },
  { day13p1, day13p2 }, { day14p1, day14p2 }, { day15p1, day15p2 },
  { day16p1, day16p2 }, { day17p1, day17p2 }, { day18p1, day18p2 },
  { day19p1, day19p2 }, { day20p1, day20p2 }, { day21p1, day21p2 },
  { day22p1, day22p2 }, { day23p1, day23p2 }, { day24p1, day24p2 },
  { day25p1, day25p2 }, { day26p1, day26p2 }, { day27p1, day27p2 },
  { day28p1, day28p2 }, { day29p1, day29p2 }, { day30p1, day30p2 },
  { day31p1, day31p2 }
};

struct AOCOpts
{
  int day = 1;
  int part = 1;
  std::string filename;
};

bool
parseArgs(int argc, char* argv[], AOCOpts& opts)
{
  int nargs = argc;
  int argind = 1;
  while (argind < nargs) {
    char* arg = argv[argind];
    // reading flag, value pairs: always look for dash first
    if (arg[0] != '-') {
      std::cerr << "Unexpected arg: " << arg << "\n";
      return false;
    }
    if (++argind == nargs) {
      std::cerr << "Expecting value following arg: " << arg << "\n";
      return false;
    }

    switch (arg[1]) {
      case 'd':
        opts.day = std::stoi(argv[argind]);
        if (opts.day <= 0 || opts.day > LAST_DAY) {
          std::cerr << "DAY must be between 1 and " << LAST_DAY << "\n";
          return false;
        }
        break;
      case 'p':
        opts.part = std::stoi(argv[argind]);
        if (opts.part != 1 && opts.part != 2) {
          std::cerr << "PART must be 1 or 2\n";
          return false;
        }
        break;
      case 'f':
        opts.filename = argv[argind];
        break;
      default:
        std::cerr << "Unexpected arg: " << arg << "\n";
        return false;
    }
    argind++;
  }
  // make it a little nicer: assume file is input/day##.txt
  if (opts.filename.empty()) {
    std::stringstream ss;
    ss <<
#ifdef WIN32
      "input\\day"
#else
      "input/day"
#endif
       << std::setfill('0') << std::setw(2) << opts.day << ".txt";
    opts.filename = ss.str();
  }
  return true;
}

int
main(int argc, char* argv[])
{
  AOCOpts opts;

  if (parseArgs(argc, argv, opts)) {
    std::cerr << "Running day " << opts.day << " part " << opts.part
              << " with input: " << opts.filename << "\n";
    std::ifstream input{ opts.filename };
    AOCResult res = (days[opts.day - 1][opts.part - 1])(input);
    std::cout << res << "\n";
    return 0;
  } else {
    return 1;
  }
}
