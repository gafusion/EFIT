#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <string.h>
#include <netdb.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>


#define PORTNUMBER 2001


void efitnotifier(void) 
{
  char localhostname[64];
  char *hostname;
  int s, len, hname;
  struct hostent *hp;
  struct sockaddr_in name;
  
  // get local host name
  hname = gethostname(localhostname, sizeof(localhostname));
  if (hname < 0) {
    perror("gethostname");
    return;
  }
  hostname = localhostname;

  // look up host's network address
  if ((hp = gethostbyname(hostname)) == NULL) {
    fprintf(stderr, "unknown host: %s. \n", hostname);
    return;
  }

  // create a socke in the INET domain
  if ((s = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
    perror("socket");
    return;
  }

  // create the address of the server
  memset(&name, 0, sizeof(struct sockaddr_in));

  name.sin_family = AF_INET;
  name.sin_port = htons(PORTNUMBER);
  memcpy(&name.sin_addr, hp->h_addr_list[0], hp->h_length);
  len = sizeof(struct sockaddr_in);

  // connect to the server
  if (connect(s, (struct sockaddr *) &name, len) < 0) {
    perror("connect - efitserver");
    return;
  }

  close(s);
  return;
}

void efitnotifier_(void)
{
void efitnotifier();
}
  
