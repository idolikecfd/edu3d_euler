/*
 * Copyright 2025. All rights reserved.
 */

#ifndef EDU3D_EULER_DEBUG_HOOKS_H
#define EDU3D_EULER_DEBUG_HOOKS_H

#ifdef __cplusplus
extern "C" {
#endif

/* Initialize debug hooks - installs signal handlers */
void edu3d_euler_debug_hooks_init(void);

/* Cleanup debug hooks - restores original signal handlers */
void edu3d_euler_debug_hooks_cleanup(void);

#ifdef __cplusplus
}
#endif

#endif /* EDU3D_EULER_DEBUG_HOOKS_H */
