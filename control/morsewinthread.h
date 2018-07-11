/**
 *
 * @file morsewinthread.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon Windows thread interface header
 *
 * @version 1.0.0
 * @author Piotr Luszczek
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 */
#ifndef CHAMELEONWINTHREAD_H
#define CHAMELEONWINTHREAD_H

#include <windows.h>

/*
typedef struct pthread_s {
  HANDLE Hth;
  unsigned IDth;
  void *(*Fth) (void *);
} pthread_t;
*/
typedef struct pthread_s {
  HANDLE hThread;
  unsigned int uThId;
} pthread_t;

typedef HANDLE pthread_mutex_t;
typedef int pthread_mutexattr_t;
typedef int pthread_attr_t;
typedef int pthread_condattr_t;

typedef struct pthread_cond_s {
  HANDLE hSem;
  HANDLE hEvt;
  CRITICAL_SECTION cs;
  int waitCount; /* waiting thread counter */
} pthread_cond_t;

typedef int pthread_attr_t;

#define PTHREAD_MUTEX_INITIALIZER ((pthread_mutex_t) -1)

#define PTHREAD_SCOPE_SYSTEM 1

#define CHAMELEON_DLLPORT
#define CHAMELEON_CDECL __cdecl

CHAMELEON_DLLPORT pthread_t CHAMELEON_CDECL pthread_self(void);
CHAMELEON_DLLPORT int CHAMELEON_CDECL pthread_mutex_init(pthread_mutex_t *mutex, const pthread_mutexattr_t * attr);
CHAMELEON_DLLPORT int CHAMELEON_CDECL pthread_mutex_destroy(pthread_mutex_t *mutex);
CHAMELEON_DLLPORT int CHAMELEON_CDECL pthread_mutex_lock(pthread_mutex_t *mutex);
CHAMELEON_DLLPORT int CHAMELEON_CDECL pthread_mutex_trylock(pthread_mutex_t *mutex);
CHAMELEON_DLLPORT int CHAMELEON_CDECL pthread_mutex_unlock(pthread_mutex_t *mutex);
CHAMELEON_DLLPORT int CHAMELEON_CDECL pthread_attr_init(pthread_attr_t *attr);
CHAMELEON_DLLPORT int CHAMELEON_CDECL pthread_attr_destroy(pthread_attr_t *attr);
CHAMELEON_DLLPORT int CHAMELEON_CDECL pthread_attr_setscope(pthread_attr_t *attr, int scope);
CHAMELEON_DLLPORT int CHAMELEON_CDECL pthread_create(pthread_t *tid, const pthread_attr_t *attr, void *(*start) (void *), void *arg);
CHAMELEON_DLLPORT int CHAMELEON_CDECL pthread_cond_init(pthread_cond_t *cond, const pthread_condattr_t *attr);
CHAMELEON_DLLPORT int CHAMELEON_CDECL pthread_cond_destroy(pthread_cond_t *cond);
CHAMELEON_DLLPORT int CHAMELEON_CDECL pthread_cond_wait(pthread_cond_t *cond, pthread_mutex_t *mutex);
CHAMELEON_DLLPORT int CHAMELEON_CDECL pthread_cond_broadcast(pthread_cond_t *cond);
CHAMELEON_DLLPORT int CHAMELEON_CDECL pthread_join(pthread_t thread, void **value_ptr);
CHAMELEON_DLLPORT int CHAMELEON_CDECL pthread_equal(pthread_t thread1, pthread_t thread2);

CHAMELEON_DLLPORT int CHAMELEON_CDECL pthread_setconcurrency (int);

CHAMELEON_DLLPORT unsigned int CHAMELEON_CDECL pthread_self_id(void);

#endif
