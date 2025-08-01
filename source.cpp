#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <random>
#include <algorithm>
#include <cmath>
#include <fstream> // For file I/O
using namespace std;

void insertion_sort(std::vector<int>& arr){
    for(int j = 1; j < arr.size(); j++){
        int key = arr[j];
        int i = j - 1;
        while(i >= 0 && arr[i] > key){
            arr[i+1] = arr[i];
            i = i - 1;
        }
        arr[i + 1] = key;
    }
}

void merge(vector<int>& arr, int p , int q, int r){
    vector<int> left(q - p + 1);
    vector<int> right(r - q);
    for(int i = 0; i < left.size(); i++){
        left[i] = arr[p + i];
    }
    for(int i = 0; i < right.size(); i++){
        right[i] = arr[q + i + 1];
    }

    int k = p, i = 0, j = 0;

    while(i < left.size() && j < right.size()){
        if(left[i] <= right[j]){
            arr[k++] = left[i++];
        }
        else{ 
            arr[k++] = right[j++];
        }
    }

    while (i < left.size()){
        arr[k++] = left[i++];
    }
    while(j < right.size()){
        arr[k++] = right[j++];
    }
}

void merge_sort(vector<int>& arr, int p, int r){
    if(p >= r){
        return;
    }
    int q = (p + r) / 2;

    merge_sort(arr, p, q);
    merge_sort(arr, q + 1, r);

    merge(arr, p, q, r);
}

int partition(vector<int>& arr, int p, int r){
    int pivot = arr[r];
    int i = p - 1;
    for(int j = p; j < r; j++){
        if(arr[j] <= pivot){
            i++;
            int temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
        }
    }
    int temp = arr[i + 1];
    arr[i + 1] = arr[r];
    arr[r] = temp;
    return i + 1;
}

int randomized_partition(vector<int>& arr, int p, int r){
    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dis(p,r);
    int j = dis(gen);

    int temp = arr[r];
    arr[r] = arr[j];
    arr[j] = temp;

    return partition(arr, p,r);
}

int randomized_select(vector<int>& arr, int p, int r, int i){
    if (p == r){
        return arr[p];
    }
    int q = randomized_partition(arr, p, r);
    int k = q - p + 1;
    if (i == k){
        return arr[q];
    }
    else if(i < k){
        return randomized_select(arr, p, q - 1, i);
    }
    else{
        return randomized_select(arr, q + 1, r, i - k);
    }
}

// Function to compute c (constant factor) for predicted runtime
double compute_c(const std::vector<double>& empirical_times, const std::vector<double>& ns, const std::vector<double>& factors){
    double numerator = 0.0;
    double denominator = 0.0;
    for(size_t i = 0; i < empirical_times.size(); i++){
        numerator += empirical_times[i] * factors[i];
        denominator += factors[i] * factors[i];
    }
    return numerator / denominator;
}

int main(){
    //seed random number generator
    srand(static_cast<unsigned int>(time(nullptr)));

    //input sizes
    vector<size_t> n_values = {10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000};

    //number of iterations
    size_t m = 5;

    std::vector<double> insertion_sort_times;
    std::vector<double> merge_sort_times;
    std::vector<double> randomized_select_times;
    std::vector<double> ns;
    std::vector<double> n_log_n; // For merge sort
    std::vector<double> n_squared; // For insertion sort

    for(size_t idx = 0; idx < n_values.size(); idx++){
        size_t n = n_values[idx];
        size_t i = floor((2*n) / 3);

        double total_time1 = 0.0, total_time2 = 0.0, total_time3 = 0.0;

        for(size_t iter = 0; iter < m; ++iter){
            std::vector<int> A(n);
            for(size_t j = 0; j < n; ++j){
                A[j] = rand();
            }

            //create copies of A
            vector<int> A_Copy1 = A;
            vector<int> A_Copy2 = A;

            //measure time for insertion sort
            auto start_time = std::chrono::high_resolution_clock::now();
            insertion_sort(A);
            auto end_time = std::chrono::high_resolution_clock::now();
            double elapsed_time = std::chrono::duration<double, std::milli>(end_time - start_time).count();
            total_time1 += elapsed_time;

            //measure time for merge sort
            auto start_time1 = std::chrono::high_resolution_clock::now();
            merge_sort(A_Copy1, 0, n - 1);
            auto end_time1 = std::chrono::high_resolution_clock::now();
            double elapsed_time1 = std::chrono::duration<double, std::milli>(end_time1 - start_time1).count();
            total_time2 += elapsed_time1;

            //measure time for randomized select
            auto start_time2 = std::chrono::high_resolution_clock::now();
            randomized_select(A_Copy2, 0, n - 1, i);
            auto end_time2 = std::chrono::high_resolution_clock::now();
            double elapsed_time2 = std::chrono::duration<double, std::milli>(end_time2 - start_time2).count();
            total_time3 += elapsed_time2;
        }

        //calculate average times
        double avg_time1 = total_time1 / m , avg_time2 = total_time2 / m , avg_time3 = total_time3 / m ;

        // Collect average times
        insertion_sort_times.push_back(avg_time1);
        merge_sort_times.push_back(avg_time2);
        randomized_select_times.push_back(avg_time3);
        ns.push_back(static_cast<double>(n));
        n_log_n.push_back(n * log(n));
        n_squared.push_back(static_cast<double>(n) * n);

        cout<< "n: " << n << endl;
        cout<< "Average time for insertion sort: " << avg_time1 << " ms" << endl;
        cout<< "Average time for merge sort: " << avg_time2 << " ms" << endl;
        cout<< "Average time for randomized select: " << avg_time3 << " ms" << endl;
        cout<< "----------------------------------------" << endl;
    }

    // Compute constants c1, c2, c3 for predicted runtimes
    double c_insertion_sort = compute_c(insertion_sort_times, ns, n_squared);
    double c_merge_sort = compute_c(merge_sort_times, ns, n_log_n);
    double c_randomized_select = compute_c(randomized_select_times, ns, ns);

    // Write data to CSV file
    ofstream outfile("runtime_data.csv");
    outfile << "n,EmpiricalRT_InsertionSort,PredictedRT_InsertionSort,EmpiricalRT_MergeSort,PredictedRT_MergeSort,EmpiricalRT_RandomizedSelect,PredictedRT_RandomizedSelect\n";
    for(size_t i = 0; i < n_values.size(); i++){
        double predicted_time_insertion_sort = c_insertion_sort * n_squared[i];
        double predicted_time_merge_sort = c_merge_sort * n_log_n[i];
        double predicted_time_randomized_select = c_randomized_select * ns[i];

        outfile << ns[i] << ",";
        outfile << insertion_sort_times[i] << "," << predicted_time_insertion_sort << ",";
        outfile << merge_sort_times[i] << "," << predicted_time_merge_sort << ",";
        outfile << randomized_select_times[i] << "," << predicted_time_randomized_select << "\n";
    }
    outfile.close();

    cout << "Data written to runtime_data.csv" << endl;

    return 0;
}
