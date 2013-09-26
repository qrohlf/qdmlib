#include <qdmlib.h>
#include <FPT.h>

#define x_size 600
#define y_size 600

int main(int argc, char **argv)
{
    FILE *f;

    if(argc == 2) {
        f = fopen(argv[1], "r");
        printf("File %s read.",argv[1]);
    } else if (argc == 1) {

    } else {
        printf("Please enter the correct amount of command line arguments.\n");
    }

    /*
     * Initialize graphics engine using sizes defined
     */
    G_init_graphics(x_size, y_size);

    /*
     * Display until a key is pressed
     */
    G_wait_key();
}
