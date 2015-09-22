#ifndef MFS_TIM_H
#define MFS_TIM_H

/**
--------------------------------------------------------------------------------
-   Module      :   mfs.cpp
-   Description :   A general purpose library for manipulating a memory area
-                   as if it were a file.
-                   mfs_ stands for memory file system.
-   Author      :   Mike Johnson - Banctec AB 03/07/96
-   Adaptations :   Tim Zaman - Picturae, Pixelprisma, 06/05/2015
--------------------------------------------------------------------------------
*/

/*

Copyright (c) 1996 Mike Johnson
Copyright (c) 1996 BancTec AB

Permission to use, copy, modify, distribute, and sell this software
for any purpose is hereby granted without fee, provided
that (i) the above copyright notices and this permission notice appear in
all copies of the software and related documentation, and (ii) the names of
Mike Johnson and BancTec may not be used in any advertising or
publicity relating to the software.

THE SOFTWARE IS PROVIDED "AS-IS" AND WITHOUT WARRANTY OF ANY KIND,
EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY
WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

Additions: Eric Bruneton, Antoine Begault, Guillaume Piolat, Tim Zaman
 
*/

typedef struct {
    char *buf = NULL;        /* Memory for open buf */
    long  buf_off;    /* File pointer for each buf */
    long  buf_size = 0;   /* Count of bytes allocated for each buf */
    int   buf_mode;   /* Mode of buffer (r, w, a) */
    bool  buf_open;
} mfs_file;

int mfs_open (void *ptr, int size, char *mode, mfs_file *fd);
int mfs_lseek (mfs_file *fd, int offset, int whence);
int mfs_read (mfs_file *fd, void *buf, int size);
int mfs_write (mfs_file *fd, void *buf, int size);
int mfs_size (mfs_file *fd);
int mfs_map (mfs_file *fd, char **addr, size_t *len);
int mfs_unmap (mfs_file *fd);
int mfs_close (mfs_file *fd);

#endif
