#include <chrono>

class Timer {
    public:
        Timer(void) {}
        
        void start(void) {
            start_time = std::chrono::steady_clock::now();
        }

        double elapsed_seconds(void) const {
            std::chrono::steady_clock::time_point now = std::chrono::steady_clock::now();
            std::chrono::duration<double> diff = now - start_time;
            return diff.count();
        }

        double elapsed_minutes(void) const {

            auto duration_mins = elapsed_seconds()/60;
            return duration_mins;
        }

        void reset(void) {
            start_time = std::chrono::steady_clock::now();
        }

        void view(std::ostream & os) {
            double time = elapsed_seconds();
            int time_m = static_cast<int>(time / 60);
            double time_s = time-(time_m*60);

            os << "time: \t" << time_m << "m" << time_s << "s";
        }

    private:
        std::chrono::steady_clock::time_point start_time;
};