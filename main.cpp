#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <cmath>
#include <chrono>

// Количество городов
const int numCities = 10;

// Количество поколений
const int numGenerations = 1000;

// Размер популяции
const int populationSize = 1000;

// Вероятность мутации
const double mutationRate = 0.02;

// Структура города
struct City {
    int x;
    int y;
};

// Расстояние между двумя городами
double distance(const City& city1, const City& city2) {
    int xDiff = city1.x - city2.x;
    int yDiff = city1.y - city2.y;
    return std::sqrt(xDiff * xDiff + yDiff * yDiff);
}

// Генерация случайного города
City generateRandomCity() {
    return {rand() % 100, rand() % 100};
}

// Генерация случайного пути (особи)
std::vector<int> generateRandomPath() {
    std::vector<int> path(numCities);
    for (int i = 0; i < numCities; ++i) {
        path[i] = i;
    }
    std::random_shuffle(path.begin(), path.end());
    return path;
}

// Расчет фитнес-функции для пути (особи)
double calculateFitness(const std::vector<int>& path, const std::vector<City>& cities) {
    double totalDistance = 0.0;
    for (int i = 0; i < numCities - 1; ++i) {
        totalDistance += distance(cities[path[i]], cities[path[i + 1]]);
    }
    // Добавляем расстояние между последним и первым городами (стартовым и конечным)
    totalDistance += distance(cities[path[numCities - 1]], cities[path[0]]);
    return totalDistance;
}


// Кроссовер двух путей (особей)
std::vector<int> crossover(const std::vector<int>& path1, const std::vector<int>& path2) {
    std::vector<int> child(numCities, -1);
    int startPos = rand() % numCities;
    int endPos = rand() % numCities;
    if (startPos > endPos) {
        std::swap(startPos, endPos);
    }
    for (int i = startPos; i <= endPos; ++i) {
        child[i] = path1[i];
    }
    int j = 0;
    for (int i = 0; i < numCities; ++i) {
        if (j == startPos) {
            j = endPos + 1;
        }
        if (std::find(child.begin(), child.end(), path2[i]) == child.end()) {
            child[j] = path2[i];
            ++j;
        }
    }
    return child;
}

// Мутация пути (особи)
void mutate(std::vector<int>& path) {
    for (int i = 0; i < numCities; ++i) {
        if (rand() / static_cast<double>(RAND_MAX) < mutationRate) {
            int j = rand() % numCities;
            std::swap(path[i], path[j]);
        }
    }
}

// Вывод лучшей длины пути
void printBestPath(const std::vector<int>& bestPath, double bestDistance, const std::vector<City>& cities) {
    std::cout << "Лучший путь: ";
    for (int i = 0; i < numCities; ++i) {
        std::cout << "(" << cities[bestPath[i]].x << "," << cities[bestPath[i]].y << ") ";
    }
    std::cout << "(" << cities[bestPath[0]].x << "," << cities[bestPath[0]].y << ")" << std::endl;
    std::cout << "Длина пути: " << bestDistance << std::endl;
}

int main() {
    // Засекаем время в начале выполнения программы

    auto start = std::chrono::high_resolution_clock::now();
    srand(time(0));
    std::ifstream inputFile("C:\\Users\\dell\\CLionProjects\\genetig\\points.txt");

    if (!inputFile) {
        std::cerr << "Не удалось открыть файл input.txt" << std::endl;
        return 1;
    }

    // Чтение координат городов из файла
    std::vector<City> cities(numCities);
    for (int i = 0; i < numCities; ++i) {
        inputFile >> cities[i].x >> cities[i].y;
    }
    inputFile.close();

    // Инициализация популяции случайными путями
    std::vector<std::vector<int>> population(populationSize);
    for (int i = 0; i < populationSize; ++i) {
        population[i] = generateRandomPath();
    }

    std::vector<double> fitness(populationSize);
    std::vector<int> offspring(numCities);


    // Основной цикл генетического алгоритма
    std::vector<int> bestPath;
    double bestDistance = std::numeric_limits<double>::max();

    bestPath.push_back(0); // 0 - индекс стартового города
    for (int g = 0; g < numGenerations; ++g) {
        // Расчет фитнес-функции для текущей популяции
        for (int i = 0; i < populationSize; ++i) {
            fitness[i] = calculateFitness(population[i], cities);
        }

        // Выбор лучшего пути в текущем поколении
        int bestIndex = std::min_element(fitness.begin(), fitness.end()) - fitness.begin();
        double currentBestDistance = fitness[bestIndex];

        if (currentBestDistance < bestDistance) {
            bestDistance = currentBestDistance;
            bestPath = population[bestIndex];
        }


        // Создание новой популяции
        for (int i = 0; i < populationSize; ++i) {
            int parentIndex1 = rand() % populationSize;
            int parentIndex2 = rand() % populationSize;
            std::vector<int>& parent1 = population[parentIndex1];
            std::vector<int>& parent2 = population[parentIndex2];

            // Кроссовер
            offspring = crossover(parent1, parent2);

            // Мутация
            mutate(offspring);

            // Замена родителей на потомка, если он лучше
            double parentFitness = fitness[parentIndex1] + fitness[parentIndex2];
            double offspringFitness = calculateFitness(offspring, cities);

            if (offspringFitness < parentFitness) {
                population[i] = offspring;
            }
            else {
                population[i] = parent1;
            }
        }

    }
    printBestPath(bestPath, bestDistance, cities);
    // Засекаем время в конце выполнения программы
    auto end = std::chrono::high_resolution_clock::now();

    // Рассчитываем время выполнения программы в миллисекундах
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    // Выводим время выполнения программы
    std::cout << "Время выполнения программы: " << duration << " мс" << std::endl;





    return 0;
}
