// ConsoleHelper.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#pragma once
#include <iostream>
#include "../Signals/rwfile.h"
#define LOGS_PATH "../Signals/Logs/"
using namespace std;
int main()
{
    File_helper fh;
    string file_name; 
    file_name= "Study2020_noizelvl_100signalSize_69998";
    file_name = LOGS_PATH + file_name;
    vector<double> vec;
    fh.read_vector_from_file(vec, file_name);
    for (auto i : vec) cout << i << endl;
}
