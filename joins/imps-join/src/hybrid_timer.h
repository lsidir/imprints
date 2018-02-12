#ifndef _HYBRID_TIMER_H_
#define _HYBRID_TIMER_H_

#include <sys/time.h>

/*
   A hybrid timer that reports both the number of CPU cycles and wall clock time (in seconds and milliseconds).
   It uses the RDTSC instruction to obtain CPU cycles,
   and uses gettimeofday() to calculate elapsed seconds.
*/

typedef unsigned long long cycle_t;

static unsigned long long rdtsc(void);

#if 0
class HybridTimer{
public:
    void Start();
    void Stop();

    cycle_t GetNumCycles() const{
    	return cycles_;
    }
    double GetSeconds() const{
    	return seconds_;
    }

    double GetMilliSeconds() const {
    	return seconds_ * 1000.0;
    }

    cycle_t GetAccuCycles() const {
    	return rdtsc() - cycles_;
    }

    double GetAccuSeconds() const {
        struct timeval endtime;
        gettimeofday(&endtime, NULL);

        return 1.0*(endtime.tv_sec - time_.tv_sec) + 1.0*(endtime.tv_usec - time_.tv_usec)/1000000;
    }

private:
    cycle_t cycles_;
    double seconds_;
    struct timeval time_;
};

inline void HybridTimer::Start(){
    cycles_ = rdtsc();
    gettimeofday(&time_, NULL);
}

inline void HybridTimer::Stop(){
    cycles_ = rdtsc() - cycles_;
    
    struct timeval endtime;
    gettimeofday(&endtime, NULL);
    seconds_ = 1.0*(endtime.tv_sec - time_.tv_sec) + 1.0*(endtime.tv_usec - time_.tv_usec)/1000000;
}
#endif

typedef struct {
    cycle_t cycles_;
    double seconds_;
    struct timeval time_;
} HybridTimer;

void timer_start(HybridTimer* timer) {
	timer->cycles_ = rdtsc();
    gettimeofday(&timer->time_, NULL);
}

void timer_stop(HybridTimer* timer) {
	timer->cycles_ = rdtsc() - timer->cycles_;

    struct timeval endtime;
    gettimeofday(&endtime, NULL);
    timer->seconds_ = 1.0*(endtime.tv_sec - timer->time_.tv_sec) + 1.0*(endtime.tv_usec - timer->time_.tv_usec)/1000000;
}

cycle_t GetNumCycles(HybridTimer* timer){
	return timer->cycles_;
}

double GetSeconds(HybridTimer* timer){
	return timer->seconds_;
}

double GetMilliSeconds(HybridTimer* timer) {
	return timer->seconds_ * 1000.0;
}

//Getting CPU cycles with RDTSC instructions
#if defined(__i386__)

static __inline__ unsigned long long rdtsc(void)
{
  unsigned long long int x;
     __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
     return x;
}
#elif defined(__x86_64__)

static __inline__ unsigned long long rdtsc(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

#elif defined(__powerpc__)

static __inline__ unsigned long long rdtsc(void)
{
  unsigned long long int result=0;
  unsigned long int upper, lower,tmp;
  __asm__ volatile(
                "0:                  \n"
                "\tmftbu   %0           \n"
                "\tmftb    %1           \n"
                "\tmftbu   %2           \n"
                "\tcmpw    %2,%0        \n"
                "\tbne     0b         \n"
                : "=r"(upper),"=r"(lower),"=r"(tmp)
                );
  result = upper;
  result = result<<32;
  result = result|lower;

  return(result);
}

#else

#error "No tick counter is available!"

#endif

#endif	//_HYBRID_TIMER_H_

