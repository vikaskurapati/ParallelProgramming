//
// Created by Dennis-Florian Herr on 13/06/2022.
//

#include <string>
#include <deque>
#include <future>
#include <functional>
#include <queue>
#include <iostream>

#include "Utility.h"

#define MEASURE_TIME true
#define NUM_THREADS 31

struct Problem
{
    Sha1Hash sha1_hash;
    int problemNum;
};

class ThreadPool
{
public:
    ThreadPool(size_t);
    template <class F, class... Args>
    auto enqueue(F &&f, Args &&...args)
        -> std::future<typename std::result_of<F(Args...)>::type>;
    ~ThreadPool();

private:
    // need to keep track of threads so we can join them
    std::vector<std::thread> workers;
    // the task queue
    std::queue<std::function<void()>> tasks;

    // synchronization
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};

// the constructor just launches some amount of workers
inline ThreadPool::ThreadPool(size_t threads)
    : stop(false)
{
    for (size_t i = 0; i < threads; ++i)
        workers.emplace_back(
            [this]
            {
                for (;;)
                {
                    std::function<void()> task;

                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->condition.wait(lock,
                                             [this]
                                             { return this->stop || !this->tasks.empty(); });
                        if (this->stop && this->tasks.empty())
                            return;
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }

                    task();
                }
            });
}

// add new work item to the pool
template <class F, class... Args>
auto ThreadPool::enqueue(F &&f, Args &&...args)
    -> std::future<typename std::result_of<F(Args...)>::type>
{
    using return_type = typename std::result_of<F(Args...)>::type;

    auto task = std::make_shared<std::packaged_task<return_type()>>(
        std::bind(std::forward<F>(f), std::forward<Args>(args)...));

    std::future<return_type> res = task->get_future();
    {
        std::unique_lock<std::mutex> lock(queue_mutex);

        // don't allow enqueueing after stopping the pool
        if (stop)
            throw std::runtime_error("enqueue on stopped ThreadPool");

        tasks.emplace([task]()
                      { (*task)(); });
    }
    condition.notify_one();
    return res;
}

// the destructor joins all threads
inline ThreadPool::~ThreadPool()
{
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for (std::thread &worker : workers)
        worker.join();
}

/*
 * TODO@Students: Implement a thread safe queue.
 * Tip: use a condition variable to make threads wait when the queue is empty and there is nothing to pop().
 * https://en.cppreference.com/w/cpp/thread/condition_variable
 */
class ProblemQueue
{
public:
    void push(Problem problem)
    {
        std::lock_guard<std::mutex> lck(mtx);
        problem_queue.push_back(problem);

        cv.notify_one();
    }

    Problem pop()
    {
        std::unique_lock<std::mutex> lck(mtx);
        while(problem_queue.empty())
        {
            cv.wait(lck);
        }
        Problem p = problem_queue.front();
        problem_queue.pop_front();
        return p;
    }

    bool empty()
    {
        return problem_queue.empty();
    }

private:
    std::deque<Problem> problem_queue;
    std::mutex mtx;
    std::condition_variable cv;
};

ProblemQueue problemQueue;

// generate numProblems sha1 hashes with leadingZerosProblem leading zero bits
// This method is intentionally compute intense so you can already start working on solving
// problems while more problems are generated
void generateProblem(int seed, int numProblems, int leadingZerosProblem)
{
    srand(seed);

    for (int i = 0; i < numProblems; i++)
    {
        std::string base = std::to_string(rand()) + std::to_string(rand());
        Sha1Hash hash = Utility::sha1(base);
        do
        {
            // we keep hashing ourself until we find the desired amount of leading zeros
            hash = Utility::sha1(hash);
        } while (Utility::count_leading_zero_bits(hash) < leadingZerosProblem);
        problemQueue.push(Problem{hash, i});
    }
}

// This method repeatedly hashes itself until the required amount of leading zero bits is found
Sha1Hash findSolutionHash(Sha1Hash hash, int leadingZerosSolution)
{
    do
    {
        // we keep hashing ourself until we find the desired amount of leading zeros
        hash = Utility::sha1(hash);
    } while (Utility::count_leading_zero_bits(hash) < leadingZerosSolution);

    return hash;
}

int main(int argc, char *argv[])
{
    int leadingZerosProblem = 8;
    int leadingZerosSolution = 11;
    int numProblems = 10000;

    // Not interesting for parallelization
    Utility::parse_input(numProblems, leadingZerosProblem, leadingZerosSolution, argc, argv);
    // Sha1Hash solutionHashes[numProblems];
    std::future<Sha1Hash> solutionHashes[numProblems];

    unsigned int seed = Utility::readInput();

#if MEASURE_TIME
    struct timespec generation_start, generation_end;
    clock_gettime(CLOCK_MONOTONIC, &generation_start);
#endif

    /*
     * TODO@Students: Generate the problem in another thread and start already working on solving the problems while the generation continues
     */
    std::thread master(&generateProblem, seed, numProblems, leadingZerosProblem);
    master.detach();

    // generateProblem(seed, numProblems, leadingZerosProblem);

#if MEASURE_TIME
    clock_gettime(CLOCK_MONOTONIC, &generation_end);
    double generation_time = (((double)generation_end.tv_sec + 1.0e-9 * generation_end.tv_nsec) - ((double)generation_start.tv_sec + 1.0e-9 * generation_start.tv_nsec));
    fprintf(stderr, "Generate Problem time:  %.7gs\n", generation_time);

    struct timespec solve_start, solve_end;
    clock_gettime(CLOCK_MONOTONIC, &solve_start);
#endif

    /*
     * TODO@Students: Create worker threads that parallelize this functionality. Add the synchronization directly to the queue
     */

    ThreadPool pool(NUM_THREADS);
    int problem_no = 1;
    while (problem_no <= numProblems)
    {
        Problem p = problemQueue.pop();
        solutionHashes[p.problemNum] = pool.enqueue(findSolutionHash, p.sha1_hash, leadingZerosSolution);
        // solutionHashes[p.problemNum] = findSolutionHash(p.sha1_hash, leadingZerosSolution);
        // std::cout << "Here" << std::endl;
        problem_no += 1;
    }

#if MEASURE_TIME
    clock_gettime(CLOCK_MONOTONIC, &solve_end);
    double solve_time = (((double)solve_end.tv_sec + 1.0e-9 * solve_end.tv_nsec) - ((double)solve_start.tv_sec + 1.0e-9 * solve_start.tv_nsec));
    fprintf(stderr, "Solve Problem time:     %.7gs\n", solve_time);
#endif

    /*
     * TODO@Students: Make sure all work has finished before calculating the solution
     * Tip: Push a special problem for each thread onto the queue that tells a thread to break and stop working
     */

    Sha1Hash solution;
    // guarantee initial solution hash data is zero
    memset(solution.data, 0, SHA1_BYTES);
    // this doesn't need parallelization. it's neglectibly fast
    for (int i = 0; i < numProblems; i++)
    {
        auto temp = solutionHashes[i].get();
        solution = Utility::sha1(solution, temp);
    }

    Utility::printHash(solution);
    printf("DONE\n");

    return 0;
}
