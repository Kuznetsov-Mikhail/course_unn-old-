// ConsoleHelper.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//
#pragma once
#include <iostream>
#include "../Signals/rwfile.h"
#define LOGS_PATH "../Signals/Logs/"
using namespace std;
int main()
{
    while (1)
    {
        File_helper fh;
        string file_name;
        cout << "File name: \n";
        cin >> file_name;
        /*file_name = "Study2020_26-04-2020-18_15_00_noizelvl_0signalSize_349990k32";*/
        file_name = LOGS_PATH + file_name;
        vector<double> vec;
        fh.read_vector_from_file(vec, file_name);
        for (auto i : vec) cout << i << endl;
    }
   
}
