/**
--------------------------------------------------------------------------------
-   Module      :   notifier.h
-   Description :   A micro wrapper for the embedded OS-X notification center.
-   Author      :   Tim Zaman, 20-NOV-2015
--------------------------------------------------------------------------------
*/

/*

Copyright (c) 2015 Tim Zaman

Permission to use, copy, modify, distribute, and sell this software
for any purpose is hereby granted without fee, provided
that (i) the above copyright notices and this permission notice appear in
all copies of the software and related documentation.

THE SOFTWARE IS PROVIDED "AS-IS" AND WITHOUT WARRANTY OF ANY KIND,
EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY
WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.
 
*/

#ifndef NOTIFIER_H
#define NOTIFIER_H

#include <string>

class Notifier {
 public:
    void notify(std::string title, std::string message){
        notifyMM(title, message);
    }
 private:
#ifdef __APPLE__ 
    // Only apple can use the notification in the .mm class
    void notifyMM(std::string title, std::string message);
#else
    // On non-apple systems, use another notification system
    void notifyMM(std::string title, std::string message){
        // @TODO Display a notification on other platforms
    }
#endif 

};

#endif
