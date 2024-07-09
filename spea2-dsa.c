#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define MAX_STRING_LENGTH 50
#define MAX_CARS 100
#define ARCHIVE_SIZE 100
#define MUTATION_RATE 0.1

typedef struct {
    char company[MAX_STRING_LENGTH];
    char model[MAX_STRING_LENGTH];
    float price;
    float mileage;
    float engine_size;
} Car;

typedef struct {
    float **data;
    int size;
} Graph;

float fitness[ARCHIVE_SIZE];  // Fitness values for each car

Car archive[ARCHIVE_SIZE];
int archive_count = 0;

Graph createGraph(int size) {
    Graph graph;
    graph.size = size;
    graph.data = (float **)malloc(size * sizeof(float *));
    for (int i = 0; i < size; i++) {
        graph.data[i] = (float *)malloc(size * sizeof(float));
    }
    return graph;
}

void destroyGraph(Graph *graph) {
    for (int i = 0; i < graph->size; i++) {
        free(graph->data[i]);
    }
    free(graph->data);
    graph->size = 0;
}

float getDistance(Graph graph, int i, int j) {
    return graph.data[i][j];
}

void setDistance(Graph *graph, int i, int j, float distance) {
    graph->data[i][j] = distance;
    graph->data[j][i] = distance; // Assuming undirected graph
}


float calculateDistance(Car car1, Car car2) {
    return sqrt(pow(car1.mileage - car2.mileage, 2) + pow(car1.engine_size - car2.engine_size, 2) + pow(car1.price - car2.price, 2));
}

int doesCar1DominateCar2(Car car1, Car car2) {
    int better_in_at_least_one = 0;
    int not_worse_in_all = 1;

    // Check each objective
    if (car1.price >car2.price) {
        better_in_at_least_one = 1;
    } else if (car1.price < car2.price) {
        not_worse_in_all = 0;
    }

    if (car1.mileage >car2.mileage) {
        better_in_at_least_one = 1;
    } else if (car1.mileage < car2.mileage) {
        not_worse_in_all = 0;
    }

    if (car1.engine_size > car2.engine_size) {
        better_in_at_least_one = 1;
    } else if (car1.engine_size < car2.engine_size) {
        not_worse_in_all = 0;
    }

    // Check dominance
    if (better_in_at_least_one && not_worse_in_all) {
        return 1;  // car1 dominates car2
    }
    return 0;  // car1 does not dominate car2
}


void updateGraphAndFitness(Car cars[], int num_cars, Graph *graph) {
    float raw_fitness[num_cars]; // Raw fitness values
    int strength[num_cars];       // Strength values

    // Step 1: Calculate dominance relationship and strength
    for (int i = 0; i < num_cars; i++) {
        strength[i] = 0; // Initialize strength for each individual
        for (int j = 0; j < num_cars; j++) {
            if (doesCar1DominateCar2(cars[i], cars[j])) {
                strength[i]++;
            }
        }
    }

    // Step 2: Calculate raw fitness
    for (int i = 0; i < num_cars; i++) {
        raw_fitness[i] = 0; // Initialize raw fitness for each individual
        for (int j = 0; j < num_cars; j++) {
            if (doesCar1DominateCar2(cars[j], cars[i])) {
                raw_fitness[i] += strength[j];
            }
        }
    }

    // Step 3: Normalize raw fitness values to obtain final fitness values
    float min_raw_fitness = INFINITY, max_raw_fitness = -INFINITY;
    for (int i = 0; i < num_cars; i++) {
        min_raw_fitness = fmin(min_raw_fitness, raw_fitness[i]);
        max_raw_fitness = fmax(max_raw_fitness, raw_fitness[i]);
    }

    // Normalize raw fitness values to obtain final fitness values (fitness[i] = (raw_fitness[i] - min_raw_fitness) / (max_raw_fitness - min_raw_fitness))
    for (int i = 0; i < num_cars; i++) {
        fitness[i] = (raw_fitness[i] - min_raw_fitness) / (max_raw_fitness - min_raw_fitness);
    }

    // Update the graph distances
    for (int i = 0; i < num_cars; i++) {
        for (int j = i + 1; j < num_cars; j++) {
            float distance = calculateDistance(cars[i], cars[j]);
            setDistance(graph, i, j, distance);
        }
    }
}

void printGraph(Graph graph) {
    printf("Graph (Adjacency Matrix):\n");
    printf("\t");
    for (int i = 0; i < graph.size; i++) {
        printf("%d\t", i);
    }
    printf("\n");

    for (int i = 0; i < graph.size; i++) {
        printf("%d\t", i);
        for (int j = 0; j < graph.size; j++) {
            printf("%.2f\t", graph.data[i][j]);
        }
        printf("\n");
    }
}

void addToArchive(Car car, Car original_population[], int num_original_cars) {
    int add_car_to_archive = 1;

    // Check dominance with the original population
    for (int i = 0; i < num_original_cars; i++) {
        if (doesCar1DominateCar2(original_population[i], car)) {
            add_car_to_archive = 0;  // Don't add car to archive
            break;
        }
        if (doesCar1DominateCar2(car, original_population[i])) {
            // Remove dominated car from original population
            for (int j = i; j < num_original_cars - 1; j++) {
                original_population[j] = original_population[j + 1];
            }
            num_original_cars--;
            i--;  // Check the current index again
        }
    }

    // Check dominance with the archive
    for (int i = 0; i < archive_count; i++) {
        if (doesCar1DominateCar2(archive[i], car)) {
            add_car_to_archive = 0;  // Don't add car to archive
            break;
        }
        if (doesCar1DominateCar2(car, archive[i])) {
            // Remove dominated car from archive
            for (int j = i; j < archive_count - 1; j++) {
                archive[j] = archive[j + 1];
            }
            archive_count--;
            i--;  // Check the current index again
        }
    }

    if (add_car_to_archive && archive_count < ARCHIVE_SIZE) {
        archive[archive_count++] = car;
    }
}


void singlePointCrossover(Car parent1, Car parent2, Car *child) {
    int crossover_point = rand() % 3;
    switch (crossover_point) {
        case 0:
            child->price = (parent1.price + parent2.price) / 2;
            child->mileage = parent1.mileage;
            child->engine_size = parent2.engine_size;
            break;
        case 1:
            child->price = parent1.price;
            child->mileage = (parent1.mileage + parent2.mileage) / 2;
            child->engine_size = parent2.engine_size;
            break;
        case 2:
            child->price = parent1.price;
            child->mileage = parent2.mileage;
            child->engine_size = (parent1.engine_size + parent2.engine_size) / 2;
            break;
    }
}

void mutation(Car *car) {
    int attribute_to_mutate = rand() % 3;
    float mutation_scale = (float)(rand() % 10 + 1) / 100; // 0.01 to 0.10
    switch (attribute_to_mutate) {
        case 0:
            car->price *= (1 + mutation_scale);
            break;
        case 1:
            car->mileage *= (1 + mutation_scale);
            break;
        case 2:
            car->engine_size *= (1 + mutation_scale);
            break;
    }
}

void binaryTournamentSelection(int *parent1, int *parent2, int num_cars) {
    if (num_cars <= 1) {
        // Handle the case when there's only one individual in the population
        *parent1 = 0;
        *parent2 = 0;
        return;
    }

    *parent1 = rand() % num_cars;
    *parent2 = rand() % num_cars;
    
    // Ensure parent1 is different from parent2
    while (*parent1 == *parent2) {
        *parent2 = rand() % num_cars;
    }
    
    // Select the parent with better fitness
    if (fitness[*parent1] > fitness[*parent2]) {
        int temp = *parent1;
        *parent1 = *parent2;
        *parent2 = temp;
    }
}


void elitism(Car cars[], float fitness[], int num_cars) {
    // Sort the archive based on fitness values
    for (int i = 0; i < archive_count - 1; i++) {
        for (int j = i + 1; j < archive_count; j++) {
            if (fitness[i] > fitness[j]) {
                // Swap fitness values
                float temp_fitness = fitness[i];
                fitness[i] = fitness[j];
                fitness[j] = temp_fitness;

                // Swap cars
                Car temp_car = archive[i];
                archive[i] = archive[j];
                archive[j] = temp_car;
            }
        }
    }
    
    // Preserve the best individuals in the archive
    int num_to_preserve = fmin(num_cars / 2, archive_count); // Ensure not to exceed the population size
    for (int i = 0; i < num_to_preserve; i++) {
        cars[i] = archive[i];
    }
   
}

void replacement(Car cars[], int num_cars) {
    // Merge the population and archive
    Car merged_population[MAX_CARS + ARCHIVE_SIZE];
    float merged_fitness[MAX_CARS + ARCHIVE_SIZE];
    int merged_count = 0;

    for (int i = 0; i < num_cars; i++) {
        merged_population[merged_count] = cars[i];
        merged_fitness[merged_count] = fitness[i];
        merged_count++;
    }

    for (int i = 0; i < archive_count; i++) {
        merged_population[merged_count] = archive[i];
        merged_fitness[merged_count] = fitness[i];
        merged_count++;
    }

    // Sort the merged population based on fitness values
    for (int i = 0; i < merged_count - 1; i++) {
        for (int j = i + 1; j < merged_count; j++) {
            if (merged_fitness[i] > merged_fitness[j]) {
                // Swap fitness values
                float temp_fitness = merged_fitness[i];
                merged_fitness[i] = merged_fitness[j];
                merged_fitness[j] = temp_fitness;

                // Swap cars
                Car temp_car = merged_population[i];
                merged_population[i] = merged_population[j];
                merged_population[j] = temp_car;
            }
        }
    }

    // Replace individuals in the population with the best ones from the merged population
    for (int i = 0; i < num_cars; i++) {
        cars[i] = merged_population[i];
        fitness[i] = merged_fitness[i];
    }
}
void findCommonSolutionAndMinDistances(Graph graph, Car archive[], int archive_count) {
    // Initialize variables to store the common solution and minimum distances
    Car common_solution;
    float min_distances[3] = {INFINITY, INFINITY, INFINITY};
    int min_distance_indices[3] = {-1, -1, -1};

    // Iterate through the archive and find the common solution with the graph
    for (int i = 0; i < archive_count; i++) {
        int found = 0;
        // Check if the current car in the archive exists in the graph
        for (int j = 0; j < graph.size; j++) {
            // If the distance in the graph is not infinity, it means there is a connection
            if (graph.data[i][j] != INFINITY) {
                // Update the common solution and set found flag
                common_solution = archive[i];
                found = 1;
                break;
            }
        }
        if (found) {
            // Calculate the distances between the common solution and other nodes in the graph
            for (int j = 0; j < graph.size; j++) {
                // If the distance in the graph is not infinity, it means there is a connection
                if (graph.data[i][j] != INFINITY) {
                    float distance = graph.data[i][j];
                    // Update the minimum distances array if the current distance is smaller than existing ones
                    for (int k = 0; k < 3; k++) {
                        if (distance < min_distances[k]) {
                            for (int l = 2; l > k; l--) {
                                min_distances[l] = min_distances[l - 1];
                                min_distance_indices[l] = min_distance_indices[l - 1];
                            }
                            min_distances[k] = distance;
                            min_distance_indices[k] = j;
                            break;
                        }
                    }
                }
            }
        }
    }

    // Print the common solution
    printf("\nCommon Solution:\n");
    printf("\nCompany\t\tModel\t\tPrice\t\tMileage\t\tEngine Size\n");
    printf("%s\t\t%s\t\t%.2f\t\t%.2f\t\t%.2f\n", 
           common_solution.company, common_solution.model, 
           common_solution.price, common_solution.mileage, 
           common_solution.engine_size);

    // Print the three nodes with minimum distances from the common solution
    printf("\nNodes with Minimum Distances from Common Solution:\n");
    printf("\nCompany\t\tModel\t\tPrice\t\tMileage\t\tEngine Size\n");
    for (int i = 0; i < 3; i++) {
    int index = min_distance_indices[i];
    printf("%s\t\t%s\t\t%.2f\t\t%.2f\t\t%.2f\n", 
           archive[index].company, archive[index].model, 
           archive[index].price, archive[index].mileage, 
         archive[index].engine_size);
}

}




void spea2(Car cars[], int num_cars, Graph *graph, Car original_population[], int num_original_cars) {
    updateGraphAndFitness(cars, num_cars, graph);  
   
    
    // Add non-dominating solutions from the original population to the archive
    for (int i = 0; i < num_original_cars; i++) {
        addToArchive(original_population[i], original_population, num_original_cars);
    }
    
    // Add non-dominating hybrid solutions to the archive
    for (int i = 0; i < num_cars; i++) {
        addToArchive(cars[i], original_population, num_original_cars);
    }
    elitism(cars, fitness, num_cars);
    replacement(cars, num_cars);

    // Selection, crossover, and mutation
    Car new_population[MAX_CARS];
    int new_population_count = 0;
    
    // Apply mating
    while (new_population_count < num_cars) {
        int parent1, parent2;
        binaryTournamentSelection(&parent1, &parent2, num_cars);

        Car child;
        singlePointCrossover(cars[parent1], cars[parent2], &child);

        if ((float)rand() / RAND_MAX < MUTATION_RATE) {
            mutation(&child);
        }

        // Add child to new population
        new_population[new_population_count++] = child;
    }

    // Replace current population with new population
    for (int i = 0; i < num_cars; i++) {
        cars[i] = new_population[i];
    }

    // Update graph and fitness with new population
    updateGraphAndFitness(cars, num_cars, graph);

    // Update archive
    for (int i = 0; i < num_cars; i++) {
        addToArchive(cars[i], original_population, num_original_cars);
    }
}

void readFromCSV(const char* filename, Car cars[], int num_cars) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        printf("Failed to open the file.\n");
        exit(1);
    }

    // Skip the header row
    char header[256];
    fgets(header, sizeof(header), file);

    // Read car data from the CSV file
    for (int i = 0; i < num_cars; i++) {
        if (fscanf(file, "%[^,],%[^,],%f,%f,%f\n",
                   cars[i].company, cars[i].model,
                   &cars[i].price, &cars[i].mileage, &cars[i].engine_size) != 5) {
            printf("Error reading data from file.\n");
            fclose(file);
            exit(1);
        }
    }
    fclose(file);
}
void removeDuplicates(Car archive[], int *archive_count) {
    for (int i = 0; i < *archive_count - 1; i++) {
        for (int j = i + 1; j < *archive_count;) { // Note: No increment of 'j' here
            if (strcmp(archive[i].company, archive[j].company) == 0 &&
                strcmp(archive[i].model, archive[j].model) == 0 &&
                archive[i].price == archive[j].price &&
                archive[i].mileage == archive[j].mileage &&
                archive[i].engine_size == archive[j].engine_size) {
                // If two cars are identical, remove one of them
                for (int k = j; k < *archive_count - 1; k++) {
                    archive[k] = archive[k + 1];
                }
                (*archive_count)--;
                // No increment of 'j' here, as the element at 'j' has been replaced
            } else {
                // Increment 'j' only if no duplicate is found
                j++;
            }
        }
    }
}



int main() {
    srand(time(NULL));

    int num_cars, num_generations;

    printf("Enter the number of cars: ");
    scanf("%d", &num_cars);

    printf("Enter the number of generations: ");
    scanf("%d", &num_generations);

    Car cars[MAX_CARS];
    Graph graph = createGraph(num_cars);
    
    readFromCSV("data.csv", cars, num_cars);

    // Store original population and its count
    Car original_population[MAX_CARS];
    int num_original_cars = num_cars;
    memcpy(original_population, cars, num_cars * sizeof(Car));

    for (int gen = 0; gen < num_generations; gen++) {
        printf("\nGeneration %d:\n", gen + 1);

        // Call spea2 with the original population and its count
        spea2(cars, num_cars, &graph, original_population, num_original_cars);
        // Remove duplicates from the archive
        removeDuplicates(archive, &archive_count);
                
        // Print the contents of the archive
        printf("Contents of the archive after generation %d:\n", gen + 1);
        printf("\nCompany\t\tModel\t\tPrice\t\tMileage\t\tEngine Size\n");
        for (int i = 0; i < archive_count; i++) {
            printf("%s\t\t%s\t\t%.2f\t\t%.2f\t\t%.2f\n", 
                   archive[i].company, archive[i].model, 
                   archive[i].price, archive[i].mileage, 
                   archive[i].engine_size);
        }
    }
    printGraph(graph);
    findCommonSolutionAndMinDistances(graph, archive, archive_count);

    destroyGraph(&graph);
    return 0;
}