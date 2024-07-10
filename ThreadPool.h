//
// Created by Hugo Minkkinen on 9.7.2024.
//

#ifndef THREADPOOL_H
#define THREADPOOL_H
#include <future>
#include <thread>
#include <queue>


class ThreadPool {
public:
    // Constructor: initializes the thread pool with a specified number of threads
    ThreadPool(size_t numThreads) : stop(false) {
        for(size_t i = 0; i < numThreads; ++i) {
            workers.emplace_back(
                [this] {
                    while(true) {
                        std::function<void()> task;
                        {
                            // Lock the queue mutex to safely access the tasks queue
                            std::unique_lock<std::mutex> lock(this->queue_mutex);
                            // Wait until there's a task or the pool is stopped
                            this->condition.wait(lock,
                                [this] { return this->stop || !this->tasks.empty(); });
                            // If the pool is stopped and there are no more tasks, exit the thread
                            if(this->stop && this->tasks.empty())
                                return;
                            // Get the next task from the queue
                            task = std::move(this->tasks.front());
                            this->tasks.pop();
                        }
                        // Execute the task
                        task();
                    }
                }
            );
        }
    }

    // Enqueue a new task to be executed by the thread pool
    template<class F, class... Args>
    auto enqueue(F&& f, Args&&... args)
        -> std::future<typename std::invoke_result<F, Args...>::type> {
        using return_type = typename std::invoke_result<F, Args...>::type;

        // Create a packaged task from the given function and arguments
        auto task = std::make_shared< std::packaged_task<return_type()> >(
                std::bind(std::forward<F>(f), std::forward<Args>(args)...)
            );

        // Get the future associated with the packaged task
        std::future<return_type> res = task->get_future();
        {
            // Lock the queue mutex to safely modify the tasks queue
            std::unique_lock<std::mutex> lock(queue_mutex);
            // If the pool is stopped, throw an exception
            if(stop)
                throw std::runtime_error("enqueue on stopped ThreadPool");
            // Add the task to the queue
            tasks.emplace([task](){ (*task)(); });
        }
        // Notify one waiting thread that a new task is available
        condition.notify_one();
        return res;
    }

    // Destructor: stops all threads and joins them
    ~ThreadPool() {
        {
            // Lock the queue mutex to safely modify the stop flag
            std::unique_lock<std::mutex> lock(queue_mutex);
            stop = true;
        }
        // Notify all waiting threads to check the stop condition
        condition.notify_all();
        // Join all worker threads
        for(std::thread &worker: workers)
            worker.join();
    }

    // Returns the number of threads in the pool
    size_t getThreadCount() const {
        return workers.size();
    }

private:
    std::vector<std::thread> workers;    // Container for worker threads
    std::queue<std::function<void()>> tasks;  // Queue of tasks to be executed
    std::mutex queue_mutex;              // Mutex to protect access to the task queue
    std::condition_variable condition;   // Condition variable for thread synchronization
    bool stop;                           // Flag to stop the thread pool
};



#endif //THREADPOOL_H
