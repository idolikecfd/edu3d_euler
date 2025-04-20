/*
 * Copyright 2025. All rights reserved.
 * 
 * Debug hooks for EDU3D_EULER - allows attaching a debugger when a crash occurs
 */

#ifdef EDU3D_EULER_ENABLE_DEBUG_HOOKS

#include "debug_hooks.h"
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <unistd.h>
#include <string.h>

/* Flag indicating whether debugging hooks are enabled */
static int debug_hooks_enabled = 0;

/* Volatile to prevent optimization - user will modify this in debugger */
static volatile int wait_for_debugger = 1;

/* Store original signal handlers */
static struct sigaction old_sigsegv;
static struct sigaction old_sigill;
static struct sigaction old_sigfpe;
static struct sigaction old_sigabrt;
static struct sigaction old_sigbus;

/* Signal handler for crashes */
static void debug_signal_handler(int sig, siginfo_t *info, void *ucontext) {
    const char *signame = "";
    
    /* Only handle if hooks are enabled */
    if (!debug_hooks_enabled)
        goto chain;
    
    /* Get signal name for display */
    switch (sig) {
        case SIGSEGV: signame = "SIGSEGV"; break;
        case SIGILL:  signame = "SIGILL"; break;
        case SIGFPE:  signame = "SIGFPE"; break;
        case SIGABRT: signame = "SIGABRT"; break;
        case SIGBUS:  signame = "SIGBUS"; break;
        default:      signame = "UNKNOWN"; break;
    }
    
    /* Print information and instructions */
    fprintf(stderr, "\n");
    fprintf(stderr, "**********************************************************\n");
    fprintf(stderr, "* EDU3D_EULER: Caught signal %d (%s)\n", sig, signame);
    fprintf(stderr, "* Process ID: %d\n", getpid());
    fprintf(stderr, "*\n");
    fprintf(stderr, "* The process is now paused so you can attach a debugger.\n");
    fprintf(stderr, "* To attach GDB: gdb -p %d\n", getpid());
    fprintf(stderr, "* Then in GDB:   set wait_for_debugger = 0\n");
    fprintf(stderr, "*                continue\n");
    fprintf(stderr, "*\n");
    fprintf(stderr, "* To terminate without debugging, send SIGKILL:\n");
    fprintf(stderr, "* kill -9 %d\n", getpid());
    fprintf(stderr, "**********************************************************\n");
    
    /* Wait for debugger to attach and set wait_for_debugger to 0 */
    while (wait_for_debugger) {
        sleep(1);
    }
    
    fprintf(stderr, "Continuing after debugger attach...\n");
    
chain:
    {
        /* Chain to original signal handler */
        struct sigaction *old_handler = NULL;
    
        switch (sig) {
            case SIGSEGV: old_handler = &old_sigsegv; break;
            case SIGILL:  old_handler = &old_sigill; break;
            case SIGFPE:  old_handler = &old_sigfpe; break;
            case SIGABRT: old_handler = &old_sigabrt; break;
            case SIGBUS:  old_handler = &old_sigbus; break;
        }
    
        if (old_handler && old_handler->sa_sigaction) {
            old_handler->sa_sigaction(sig, info, ucontext);
        } else if (old_handler && old_handler->sa_handler) {
            if (old_handler->sa_handler == SIG_DFL) {
                /* Restore default handler and re-raise signal */
                signal(sig, SIG_DFL);
                raise(sig);
            } else if (old_handler->sa_handler != SIG_IGN) {
                old_handler->sa_handler(sig);
            }
        }
    }   
}

/* Initialize debug hooks */
void edu3d_euler_debug_hooks_init() {
    struct sigaction sa;
    
    debug_hooks_enabled = 1;
    
    /* Set up signal handler */
    memset(&sa, 0, sizeof(sa));
    sa.sa_sigaction = debug_signal_handler;
    sa.sa_flags = SA_SIGINFO;
    
    /* Install signal handlers and save old ones */
    sigaction(SIGSEGV, &sa, &old_sigsegv);
    sigaction(SIGILL, &sa, &old_sigill);
    sigaction(SIGFPE, &sa, &old_sigfpe);
    sigaction(SIGABRT, &sa, &old_sigabrt);
    sigaction(SIGBUS, &sa, &old_sigbus);
    
    fprintf(stderr, "EDU3D_EULER debug hooks initialized (PID: %d)\n", getpid());
}

/* Cleanup debug hooks */
void edu3d_euler_debug_hooks_cleanup() {
    if (!debug_hooks_enabled)
        return;
        
    /* Restore original signal handlers */
    sigaction(SIGSEGV, &old_sigsegv, NULL);
    sigaction(SIGILL, &old_sigill, NULL);
    sigaction(SIGFPE, &old_sigfpe, NULL);
    sigaction(SIGABRT, &old_sigabrt, NULL);
    sigaction(SIGBUS, &old_sigbus, NULL);
    
    debug_hooks_enabled = 0;
}

#else // EDU3D_EULER_ENABLE_DEBUG_HOOKS

void edu3d_euler_debug_hooks_init() { }

void edu3d_euler_debug_hooks_cleanup() { }

#endif // EDU3D_EULER_ENABLE_DEBUG_HOOKS
