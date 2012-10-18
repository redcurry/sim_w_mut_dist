// Numeric simulation of evolving organisms where
// the distribution of fitness effects is given by a file

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>

using namespace std;


#define N 1000    // Population size
#define U 0.1     // Mutation rate (per genome per generation)

// Relative to the ancestor, fitness of the organism
// in which the mutational distribution is based on
#define MUT_DIST_NEUTRAL_FIT 0.0765346

#define INITIAL_FITNESS 0.5    // Initial fitness of entire population
#define UPDATES 100000         // Number of replications to perform

#define REVERTANT_MIN_FIT 0.999999   // Minimum fitness of revertant


struct organism
{
  double fitness;
  int generation;
  int mutations;
  bool revertant;
};


// File functions
template<typename T> vector<T> read_values(string const & path);

// Population functions
vector<organism> create_population();
void print_pop_fitnesses(vector<organism> const & pop);
void print_header();
void print_info(vector<organism> const & pop, int rep_num);
double mean_generation(vector<organism> const & pop);
double mean_fitness(vector<organism> const & pop);
double mean_mutations(vector<organism> const & pop);
double max_fitness(vector<organism> const & pop);
int count_revertants_by_inheritance(vector<organism> const & pop);
int count_revertants_by_fitness(vector<organism> const & pop);
vector<double> get_fitnesses(vector<organism> const & pop);
vector<int> get_generations(vector<organism> const & pop);
vector<int> get_mutations(vector<organism> const & pop);

// Reproduction functions
void replicate_next_organism(vector<organism> & pop,
  vector<double> const & mut_dist);
int pick_next_position(vector<organism> const & pop);
int pick_random_position();
bool should_mutate();
double calc_new_fitness(double fitness, vector<double> const & mut_dist);
bool gained_reversion(organism org);

// Random functions
double rand_0_to_1();
template<typename T> T random_choice(vector<T> list);
int roulette_choice(vector<double> const & weights);

// General list functions
template<typename T> void print(vector<T> list);
template<typename T> T max(vector<T> const & list);
template<typename T> T sum(vector<T> const & list);
template<typename T> double mean(vector<T> const & list);


int main(int argc, char* argv[])
{
  if (argc < 4)
  {
    cout << "Arguments: random_seed mut_dist_path rep_num" << endl;
    return EXIT_FAILURE;
  }

  // Initialize command-line arguments
  int random_seed = atoi(argv[1]);
  string mut_dist_path = string(argv[2]);
  int rep_num = atoi(argv[3]);

  srand(random_seed);

  // Load the mutational distribution
  vector<double> const & mut_dist = read_values<double>(mut_dist_path);

  vector<organism> pop = create_population();
  
  // Print out header only for the first replicate
  if (rep_num == 1)
    print_header();

  for(int update = 0; update < UPDATES; update++)
  {
    replicate_next_organism(pop, mut_dist);

    if (update % 1000 == 0)
      print_info(pop, rep_num);
  }

  // print_pop_fitnesses(pop);

  return 0;
}


template<typename T>
vector<T> read_values(string const & path)
{
  vector<T> values;
  ifstream in(path.c_str());

  T value;
  while (in >> value)
    values.push_back(value);

  in.close();
  return values;
}


vector<organism> create_population()
{
  return vector<organism>(N, (organism){INITIAL_FITNESS, 0, 0, false});
}


void print_header()
{
  cout << "Replicate MeanGeneration MeanFitness MaxFitness "
    << "MeanMutations N_RevertantsByInheritance N_RevertantsByFitness\n";
}


void replicate_next_organism(vector<organism> & pop,
  vector<double> const & mut_dist)
{
  int parent_pos = pick_next_position(pop);
  int child_pos = pick_random_position();

  organism parent = pop[parent_pos];
  pop[child_pos] = parent;

  if (should_mutate())
  {
    pop[child_pos].fitness = calc_new_fitness(parent.fitness, mut_dist);
    pop[child_pos].mutations = parent.mutations + 1;
    if (gained_reversion(pop[child_pos]))
      pop[child_pos].revertant = true;
  }

  pop[child_pos].generation = parent.generation + 1;
}


int pick_next_position(vector<organism> const & pop)
{
  return roulette_choice(get_fitnesses(pop));
}


int roulette_choice(vector<double> const & weights)
{
  double frac = rand_0_to_1() * sum(weights);
  double cum_sum = weights[0];

  int i = 1;
  while (cum_sum < frac)
    cum_sum += weights[i++];

  return i - 1;
}


double rand_0_to_1()
{
  return (double)rand() / RAND_MAX;
}


template<typename T>
T sum(vector<T> const & list)
{
  T sum = 0;
  for(int i = 0; i < list.size(); i++)
    sum += list[i];
  return sum;
}


vector<double> get_fitnesses(vector<organism> const & pop)
{
  vector<double> fitnesses;
  for (int i = 0; i < pop.size(); i++)
    fitnesses.push_back(pop[i].fitness);
  return fitnesses;
}


int pick_random_position()
{
  return rand() % N;
}


bool should_mutate()
{
  return rand_0_to_1() < U;
}


double calc_new_fitness(double fitness, vector<double> const & mut_dist)
{
  double mut_effect = random_choice(mut_dist);

  if (mut_effect < MUT_DIST_NEUTRAL_FIT)
    return fitness * mut_effect / MUT_DIST_NEUTRAL_FIT;
  else
    return fitness + (1 - fitness) * (mut_effect - MUT_DIST_NEUTRAL_FIT) /
      (1 - MUT_DIST_NEUTRAL_FIT);
}


template<typename T>
T random_choice(vector<T> list)
{
  return list[rand() % list.size()];
}


bool gained_reversion(organism org)
{
  return !org.revertant && org.fitness > REVERTANT_MIN_FIT;
}


void print_info(vector<organism> const & pop, int rep_num)
{
  cout << rep_num << " " << mean_generation(pop) << " "
    << mean_fitness(pop) << " "
    << max_fitness(pop) << " "
    << mean_mutations(pop) << " "
    << count_revertants_by_inheritance(pop) << " "
    << count_revertants_by_fitness(pop) << endl;
}


double mean_generation(vector<organism> const & pop)
{
  return mean(get_generations(pop));
}


template<typename T>
double mean(vector<T> const & list)
{
  return (double)sum(list) / list.size();
}


vector<int> get_generations(vector<organism> const & pop)
{
  vector<int> generations;
  for (int i = 0; i < pop.size(); i++)
    generations.push_back(pop[i].generation);
  return generations;
}


double mean_fitness(vector<organism> const & pop)
{
  return mean(get_fitnesses(pop));
}


double max_fitness(vector<organism> const & pop)
{
  return max(get_fitnesses(pop));
}


template<typename T>
T max(vector<T> const & list)
{
  T max = 0;
  for (int i = 0; i < list.size(); i++)
    if (list[i] > max)
      max = list[i];
  return max;
}


double mean_mutations(vector<organism> const & pop)
{
  return mean(get_mutations(pop));
}


vector<int> get_mutations(vector<organism> const & pop)
{
  vector<int> mutations;
  for (int i = 0; i < pop.size(); i++)
    mutations.push_back(pop[i].mutations);
  return mutations;
}


int count_revertants_by_inheritance(vector<organism> const & pop)
{
  int sum = 0;
  for (int i = 0; i < pop.size(); i++)
    if (pop[i].revertant)
      sum++;
  return sum;
}


int count_revertants_by_fitness(vector<organism> const & pop)
{
  int sum = 0;
  for (int i = 0; i < pop.size(); i++)
    if (pop[i].fitness > REVERTANT_MIN_FIT)
      sum++;
  return sum;
}


void print_pop_fitnesses(vector<organism> const & pop)
{
  print(get_fitnesses(pop));
}


template<typename T>
void print(vector<T> list)
{
  for (int i = 0; i < list.size(); i++)
    cout << list[i] << endl;
}
