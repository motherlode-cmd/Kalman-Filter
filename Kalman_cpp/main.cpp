#include  "./model/Navigation.hpp"
#include  "./model/Detection.hpp"
#include "./model/Visualize.hpp"
#include <fstream>
#include "unistd.h"
#include <iostream>
#include <thread>
#include <atomic>
#include <condition_variable>
Navigation CommandCenter;

void Init() {
    for (int i = 0; i < 12; i++) {
        my_Vector vec(3);
        for(int i = 0; i < 3; i++) {
            vec.set(i, (rand() % 100000) | 0);
        }
        CommandCenter.createReceiver(vec);
    }
    my_Vector aircraftVec(3);
    for(int i = 0; i < 3; i++) {
        aircraftVec.set(i, rand() % 100000);
    }
    CommandCenter.createAircraft(aircraftVec);
    CommandCenter.startAircraft();
}

std::atomic<bool> isStop(false);

// Функция, выполняющаяся в отдельном потоке
void detectionThread() {
    std::ofstream real;
    std::ofstream mls_result;
    std::ofstream kalman_result;
    std::ofstream error;
    real.open("plot/real.txt");
    mls_result.open("plot/mls.txt");
    kalman_result.open("plot/filter.txt");
    error.open("plot/error.txt");
    while (!isStop.load()) {
        detection(true, CommandCenter, real, mls_result, kalman_result, error);

        std::this_thread::sleep_for(std::chrono::milliseconds(100)); // Задержка между итерациями цикла
    }
    real.close();
    mls_result.close();
    kalman_result.close();
    error.close();
}

int main(){
    std::ofstream f;
    f.open("time.txt");
    std::chrono::milliseconds currentTime = std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch()
    );
    int64_t currentTimeInMilliseconds = currentTime.count();
    f<<currentTimeInMilliseconds/1000<<std::endl;
    f.close();
    Init();

    std::thread detection(detectionThread);

    // Ожидание нажатия клавиши Enter
    std::cout << "Press Enter to stop..." << std::endl;
    std::cin.ignore();

    // Установка флага остановки
    isStop.store(true);

    // Ожидание завершения потока
    detection.join();
    CommandCenter.stopAircraft();
    // Дополнительная логика после остановки обработки сигналов
    // ...
    return 0;
}

