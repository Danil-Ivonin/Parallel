#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Определение структуры Person
typedef struct {
    char* name;
    int age;
} Person;

// Функция для создания нового Person
Person* newPerson(const char* name, int age) {
    Person* p = (Person*)malloc(sizeof(Person));
    p->name = strdup(name);
    p->age = age;
    return p;
}

// Функция для освобождения памяти, выделенной для Person
void freePerson(Person* p) {
    free(p->name);
    free(p);
}

// Функция для добавления элемента в массив
Person** append(Person** array, int* size, Person* p) {
    // Увеличиваем размер массива на 1
    *size += 1;
    array = (Person**)realloc(array, (*size) * sizeof(Person*));
    array[*size - 1] = p;
    return array;
}

int main() {
    // Инициализация массива и его размера
    Person** arr = NULL;
    int size = 0;

    // Добавление объектов в массив
    arr = append(arr, &size, newPerson("jon", 20));
    arr = append(arr, &size, newPerson("den", 22));
    arr = append(arr, &size, newPerson("kate", 10));

    // Печать имен всех объектов в массиве
    for (int i = 0; i < size; i++) {
        printf("%s\\n", arr[i]->name);
    }

    // Освобождение памяти
    for (int i = 0; i < size; i++) {
        freePerson(arr[i]);
    }
    free(arr);

    return 0;
}