/**
--------------------------------------------------------------------------------
-   Module      :   notifier.mm
-   Description :   A micro wrapper for the embedded OS-X notification center.
-   Author      :   Tim Zaman, 20-NOV-2015
--------------------------------------------------------------------------------
*/

#include "notifier.h"
#include <NSUserNotification.h>

void Notifier::notifyMM(std::string title, std::string message) {
	NSUserNotification *userNotification = [[[NSUserNotification alloc] init] autorelease];
	userNotification.title =  [NSString stringWithCString:title.c_str()  encoding:[NSString defaultCStringEncoding]];
	userNotification.informativeText =  [NSString stringWithCString:message.c_str()  encoding:[NSString defaultCStringEncoding]];
	NSUserNotificationCenter* center = [NSUserNotificationCenter defaultUserNotificationCenter];
	[center deliverNotification:userNotification];
}