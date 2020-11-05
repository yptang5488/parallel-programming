#include <iostream>
#include <random>
#include <stdlib.h>
#include <pthread.h>

using namespace std;

pthread_mutex_t MUTEX = PTHREAD_MUTEX_INITIALIZER;
long long IN_CIRCLE;

void *child(void *arg)
{
    long long toss_num = (long long)arg;
    long long in_circle = 0;

    random_device rd;
    /* 梅森旋轉演算法 */
    mt19937 generator(rd());
    uniform_real_distribution<float> unif(0.0, 1.0);

    for (long long toss = 0; toss < toss_num; toss++)
    {
        float x = unif(generator);
        float y = unif(generator);
        if (x * x + y * y <= 1) in_circle++;
    }

    pthread_mutex_lock(&MUTEX);
    IN_CIRCLE += in_circle;
    pthread_mutex_unlock(&MUTEX);
    pthread_exit(NULL);
}

void split(int core_num, long long remain, long long *info)
{
    for (int i = 0; i < core_num; i++)
    {
        long long pack = remain / (core_num - i);
        *(info + i) = pack;
        remain -= pack;
    }
    return;
}

int main(int argc, char *argv[])
{
    // string to long and long long
    long core_num = strtol(argv[1], NULL, 10);
    long long toss_num = strtoll(argv[2], NULL, 10);

    // build pthread
    pthread_t *cores = new pthread_t[core_num];
    long long *child_toss = new long long[core_num];

    split(core_num, toss_num, child_toss);
    // long long remain = toss_num;
    // int cut = core_num; 
    // for(int i = 0 ; i < core_num ; i++){
    //     long long pack = remain / cut;
    //     *(child_toss + i) = pack;
    //     remain -= pack;
    //     cut--;
    // }

    for (int i = 0; i < core_num; i++)
    {
        // cout << "child toss " << child_toss[i] << endl;
        pthread_create(&cores[i], NULL, child, (void *)child_toss[i]);
    }
    for (int i = 0; i < core_num; i++)
    {
        pthread_join(cores[i], NULL);
    }
    cout << 4 * IN_CIRCLE / (double)toss_num << endl;
    return 0;
} 