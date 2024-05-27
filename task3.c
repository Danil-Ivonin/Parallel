#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <stdbool.h>
#include <ctype.h>

#define TABLE_SIZE 5000

typedef struct Node {
    char* word;
    int count;
    struct Node* next;
} Node;

unsigned int hash(const char* word) 
{
    unsigned int value = 0;
    for (; *word; ++word) 
    {
        value = value * 512 + *word;
    }
    return value % TABLE_SIZE;
}

Node* create_node(const char* word) 
{
    Node* new_node = (Node*)malloc(sizeof(Node));
    new_node->word = strdup(word);
    new_node->count = 1;
    new_node->next = NULL;
    return new_node;
}

void add_word(Node** hash_table, const char* word) 
{
    unsigned int index = hash(word);
    Node* list = hash_table[index];
    for (Node* node = list; node != NULL; node = node->next) 
    {
        if (strcmp(node->word, word) == 0) 
        {
            node->count++;
            return;
        }
    }
    Node* new_node = create_node(word);
    new_node->next = list;
    hash_table[index] = new_node;
}

// Функция для сравнения двух узлов по частоте слов
int compare_nodes(const void* a, const void* b) {
    Node* node1 = *(Node**)a;
    Node* node2 = *(Node**)b;
    return node2->count - node1->count; // Сортировка по убыванию частоты
}

// Функция для сортировки hash_table
Node** sort_hash_table(Node** hash_table) {
    // Сначала подсчитаем количество слов
    int word_count = 0;
    for (int i = 0; i < TABLE_SIZE; ++i) {
        for (Node* node = hash_table[i]; node != NULL; node = node->next) {
            word_count++;
        }
    }

    // Создаем массив указателей на узлы
    Node** nodes_array = malloc(word_count * sizeof(Node*));
    int index = 0;
    for (int i = 0; i < TABLE_SIZE; ++i) {
        for (Node* node = hash_table[i]; node != NULL; node = node->next) {
            nodes_array[index++] = node;
        }
    }

    // Сортируем массив узлов
    qsort(nodes_array, word_count, sizeof(Node*), compare_nodes);
    return nodes_array;
}

void free_table(Node** hash_table) {
    for (int i = 0; i < TABLE_SIZE; ++i) {
        Node* node = hash_table[i];
        while (node != NULL) {
            Node* temp = node;
            node = node->next;
            free(temp->word);
            free(temp);
        }
    }
}

long int findSize(const char* filename) 
{ 
    // opening the file in read mode 
    FILE* fp = fopen(filename, "r"); 
  
    // checking if the file exist or not 
    if (fp == NULL) { 
        printf("File Not Found!\n"); 
        return -1; 
    } 
  
    fseek(fp, 0L, SEEK_END); 
  
    // calculating the size of the file 
    long int res = ftell(fp); 
  
    // closing the file 
    fclose(fp); 
  
    return res; 
} 

int readtext(long int text_len, char* out_text , const char* filename)
{

    FILE *fp;
    fp = fopen (filename, "r");
    if (fp == NULL)
    {
        printf("Error: could not open file %s", filename);
        return 1;
    }

    for(long int i = 0; i < text_len; ++i)
    {
        out_text[i] = fgetc(fp);
    }


    fclose (fp); 
    return 1; 
}

int main(int argc, char** argv) 
{
    int rank, size;
    double tstart, tend;char* text_name = "vishnevyi-sad";
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //char* text_name = "vishnevyi-sad";
    char* text_name = "shinel";
    char filename[100];
    sprintf(filename, "%s.txt", text_name); 
    long int text_len = findSize(filename);

    long int chunk = ceil((double)text_len / size);
    long int start = rank * chunk;
    long int end = (rank+1) * chunk - 1;
    if (end > text_len)
        end = text_len - 1;
    
    char* text = (char*)malloc(text_len * sizeof(char));
    if (rank == 0)
    {
        readtext(text_len, text, filename);
    }

    char* text_part = (char*)malloc(chunk * sizeof(char));
    MPI_Scatter(text, chunk, MPI_CHAR, text_part, chunk, MPI_CHAR, 0, MPI_COMM_WORLD);    
    tstart = MPI_Wtime();

    char text_part_copy[chunk];
    strcpy(text_part_copy, text_part);

    Node* hash_table[TABLE_SIZE] = {NULL};
    char* word = strtok(text_part_copy, " ,.!?():;\n"); // Разделители слов

    while (word != NULL) 
    {
        add_word(hash_table, word);
        word = strtok(NULL, " ,.!?():;\n"); // Получаем следующее слово
    }
    Node** sorted_table = sort_hash_table(hash_table);

    sprintf(filename, "%s_%d_words.txt", text_name, rank);
    FILE *pf;
    pf = fopen(filename, "w+");
    if (pf == NULL)
        return 0;

    for(size_t i = 0; i < 10; ++i)
    {
        fprintf(pf, "%s:%d\n", sorted_table[i]->word, sorted_table[i]->count);
    }

    fclose (pf); 
    free_table(hash_table);
    MPI_Barrier(MPI_COMM_WORLD);

    double mean_sentence_len = 0.0;
    int sent_counter = 0;

    char* sentence = strtok(text_part, ".!?;");
    while (sentence != NULL) 
    {
        if (strlen(sentence) > 1)
        {
            sent_counter += 1;
        }
        sentence = strtok(NULL, ".!?;"); // Получаем следующее слово
    }

    mean_sentence_len = (chunk - sent_counter) / sent_counter;

    double* mean_sentence_lens = (double*)malloc(size * sizeof(double));
    MPI_Gather(&mean_sentence_len, 1, MPI_DOUBLE, mean_sentence_lens, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        double mean_len_sum = 0.0;
        for (int i = 0; i < size; i++)
            mean_len_sum += mean_sentence_lens[i];
        mean_sentence_len = mean_len_sum / size;

        sprintf(filename, "%s_data.txt", text_name);
        FILE *pf;
        pf = fopen(filename, "w+");
        if (pf == NULL)
            return 0;
        fprintf(pf, "%d\n", size);
        fprintf(pf, "%f", mean_sentence_len);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    tend = MPI_Wtime();
    if (rank == 0)
    {
        double time_taken_parallel = tend - tstart;
        printf("Time spend = %f\n", time_taken_parallel);
    }
    MPI_Finalize();
    return 0;
}