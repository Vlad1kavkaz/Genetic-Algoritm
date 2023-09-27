#include <iostream>
#include <random>
#include <vector>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <map>
#include <chrono>
#include <iomanip>
#include <numeric>

const int POPULATOIN_SIZE = 20;
const int NUM_GENERATIONS = 5000;
const double MUTATION_RATE = 0.08;
const int NUM_CITIES = 11;
const double BEST_DISTANCE = 20.8438;

struct Point
{
    int xCoord;
    int yCoord;

    Point(int x, int y) : xCoord(x), yCoord(y) {}
    Point() : xCoord(0), yCoord(0) {}
};

bool operator == (const Point point1, Point point2)
{
    bool f = false;
    if (point1.xCoord == point2.xCoord && point1.yCoord == point2.yCoord)
    {
        f = true;
    }
    return f;
}

std::vector<Point> read_points_from_file(const std::string& filename)
{

    std::vector<Point> points (NUM_CITIES);
    std::ifstream in(filename);
    if (in.is_open())
    {
        for (int i = 0; i < NUM_CITIES; ++i)
        {
            int _x, _y;
            in >> _x >>_y;
            Point P(_x, _y);
            points[i] = P;
        }

    }
    in.close();
    return points;
}

double distance (Point point1, Point point2)
{
    double distance;
    distance = sqrt(pow(point2.xCoord - point1.xCoord, 2) + pow(point2.yCoord - point1.yCoord, 2));
    return distance;
}

std::vector<std::vector<int>> create_population() {
    std::vector<std::vector<int>> population(POPULATOIN_SIZE);
    std::vector<int> indices(NUM_CITIES);
    std::iota(indices.begin(), indices.end(), 0);

    auto currentTime = std::chrono::system_clock::now().time_since_epoch();
    unsigned seed = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime).count();
    std::mt19937 generator(seed);

    for (int i = 0; i < POPULATOIN_SIZE; ++i) {
        std::vector<int> individual = indices;

        std::shuffle(individual.begin(), individual.end(), generator);

        population[i] = individual;
    }

    return population;
}

double path_length(const std::vector<int>& points_order, const std::vector<Point>& points)
{
    double total_distance = 0.0;

    for (int index = 0; index < NUM_CITIES - 1; ++index)
    {
        int index1 = points_order[index];
        int index2 = points_order[index + 1];
        if (index1 < 0 || index1 >= NUM_CITIES || index2 < 0 || index2 >= NUM_CITIES) {
            std::cout << "Invalid result: index out of range." << std::endl;
            return total_distance;
        }
        Point point1 = points[index1];
        Point point2 = points[index2];
        total_distance += distance(point1, point2);
    }
    int first_index = points_order.front();
    int last_index = points_order.back();
    if (first_index >= 0 && first_index < NUM_CITIES && last_index >= 0 && last_index < NUM_CITIES) {
        total_distance += distance(points[last_index], points[first_index]);
    }
    return total_distance;
}

std::vector<int> chooseIndices(const std::vector<double>& probabilities, int num_indices) {
    std::vector<int> selected_indices(num_indices);

    auto currentTime = std::chrono::system_clock::now().time_since_epoch();
    unsigned seed = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime).count();
    std::mt19937 gen(seed);
    std::discrete_distribution<> dis(probabilities.begin(), probabilities.end());

    for (int i = 0; i < num_indices; ++i) {
        selected_indices[i] = dis(gen);
    }

    return selected_indices;
}

std::vector<std::vector<int>> selection(const std::vector<std::vector<int>>& population, const std::vector<Point>& points) {
    std::vector<double> fitness_scores(population.size());
    double total_fitness = 0.0;

    for (int i = 0; i < population.size(); ++i) {
        fitness_scores[i] = 1 / path_length(population[i], points);
        total_fitness += fitness_scores[i];
    }

    std::vector<double> probabilities(population.size(), 0.0);

    double sum = 0.0;
    for (int i = 0; i < population.size(); ++i) {
        probabilities[i] = sum += fitness_scores[i] / total_fitness;
    }

    std::vector<int> selected_indices = chooseIndices(probabilities, 2);
    std::vector<std::vector<int>> result_population(selected_indices.size());

    for (int i = 0; i < selected_indices.size(); ++i) {
        result_population[i] = population[selected_indices[i]];
    }

    return result_population;
}

std::vector<int> crossover(const std::vector<int>& parent1, const std::vector<int>& parent2) {
    auto currentTime = std::chrono::system_clock::now().time_since_epoch();
    unsigned seed = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime).count();
    std::mt19937 gen(seed);
    std::uniform_int_distribution<size_t> dis_size(0, parent1.size() - 1);
    std::uniform_int_distribution<size_t> dis_i(0, parent1.size() - 1);

    size_t i = dis_i(gen);
    size_t j = dis_size(gen);

    if (i > j) {
        std::swap(i, j);
    }

    std::vector<int> child(parent1.size(), -1);
    std::vector<bool> used(child.size(), false);

    for (size_t index = i; index < j; ++index) {
        child[index] = parent1[index];
        used[parent1[index]] = true;
    }

    size_t child_index = 0;

    for (size_t index = 0; index < parent2.size(); ++index) {
        if (child[index] == -1) {
            while (used[parent2[child_index]]) {
                ++child_index;
                child_index %= parent2.size();
            }
            child[index] = parent2[child_index];

            used[parent2[child_index]] = true;
            ++child_index;
            child_index %= parent2.size();
        }
    }

    return child;
}

std::vector<int> mutation(std::vector<int> child) {
    auto currentTime = std::chrono::system_clock::now().time_since_epoch();
    unsigned seed = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime).count();
    std::mt19937 gen(seed);
    std::uniform_int_distribution<size_t> dis(0, child.size() - 1);

    size_t mut_count = child.size() / 2;
    for (size_t i = 0; i < mut_count; ++i) {
        size_t index1 = dis(gen);
        size_t index2 = dis(gen);
        std::swap(child[index1], child[index2]);
    }

    return child;
}


std::vector<int> opt_2(std::vector<int> path, const std::vector<Point>& points)
{

    bool improved = true;
    while (improved) {
        improved = false;
        for (int i = 1; i < NUM_CITIES - 1; ++i) {
            for (int j = i + 1; j < NUM_CITIES - 1; ++j) {
                double distance1 = distance(points[path[i - 1]], points[path[i]]);
                double distance2 = distance(points[path[j]], points[path[j + 1]]);
                double distance3 = distance(points[path[i - 1]], points[path[j]]);
                double distance4 = distance(points[path[i]], points[path[j + 1]]);
                double current_distance = distance1 + distance2;
                double new_distance = distance3 + distance4;

                if (new_distance < current_distance) {
                    std::reverse(path.begin() + i, path.begin() + j + 1);
                    improved = true;
                }
            }
        }
    }

    return path;
}

std::map<std::vector<int>, double> genetic_algorithm(const std::vector<Point>& points)
{
    std::vector<std::vector<int>> population = create_population();
    double best_distance = std::numeric_limits<double>::infinity();
    std::vector<int> best_path;

    auto currentTime = std::chrono::system_clock::now().time_since_epoch();
    unsigned seed = std::chrono::duration_cast<std::chrono::milliseconds>(currentTime).count();
    std::mt19937 engine(seed);
    std::uniform_real_distribution<double> distribution(0.0, 1.0);

    for (int i = 0; i < NUM_GENERATIONS; ++i) {
        std::vector<std::vector<int>> selected_population = selection(population, points);
        std::vector<std::vector<int>> new_population;
        for (size_t j = 0; j < selected_population.size() - 1; j += 2) {
            std::vector<int> child1 = selected_population[j];
            std::vector<int> child2 = selected_population[j + 1];

            for (const auto& individual : population) {
                new_population.emplace_back(individual);
            }
            if (distribution(engine) < MUTATION_RATE) {
                child1 = mutation(child1);
            }
            if (distribution(engine) < MUTATION_RATE) {
                child2 = mutation(child2);
            }

            auto child1Iterator = std::find(new_population.begin(), new_population.end(), child1);
            if (child1Iterator != new_population.end()) {
                new_population.erase(child1Iterator);
            }

            std::vector<int> child = crossover(child1, child2);
            new_population.emplace_back(child);
        }

        population = new_population;

        std::vector<int> current_best_path = population[0];
        double min_path_length = path_length(current_best_path, points);

        for (const auto& individual : population) {
            double length = path_length(individual, points);
            if (length < min_path_length) {
                min_path_length = length;
                current_best_path = individual;
            }
        }
        double current_best_distance = path_length(current_best_path, points);
        if (current_best_distance < best_distance) {
            best_path = current_best_path;
            best_distance = current_best_distance;
        }
        if (best_distance <= BEST_DISTANCE)
        {
            break;
        }
    }

    best_path.push_back(best_path[0]);
    best_path = opt_2(best_path, points);

    std::map<std::vector<int>, double> result;
    result.insert(std::make_pair(best_path, best_distance));
    return result;
}

void print_result(std::map<std::vector<int>, double> result, std::vector<Point> points)
{
    std::vector<int> path = result.begin()->first;
    std::ofstream out ("C:\\Users\\28218\\CLionProjects\\genetic\\output.txt");
    std::cout << "Best path: ";
    for (int i = 0; i < path.size() - 1; ++i)
    {
        std::cout << "(" << points[path[i]].xCoord << ";" << points[path[i]].yCoord << ") -> ";
        out << points[path[i]].xCoord << " " << points[path[i]].yCoord << "\n";
    }
    std::cout << "(" << points[path.back()].xCoord << ";" << points[path.back()].yCoord << ")" << std::endl;
    out<< points[path.back()].xCoord << " " << points[path.back()].yCoord;
    out.close();
    std::cout << "Best distances: ";
    std::cout << result.begin()->second;
}

int main()
{
    std::vector<Point> points = read_points_from_file("C:\\Users\\28218\\CLionProjects\\genetic\\points.txt");
    auto start = std::chrono::high_resolution_clock::now();
    std::map<std::vector<int>, double> result = genetic_algorithm(points);
    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    double duration_seconds = duration / 1000.0;
    std::cout << "Program running time: " << std::fixed << std::setprecision(5) << duration_seconds << " seconds" << std::endl;
    print_result(result, points);

    return 0;
}
