# crow-example

# An example of NMTF using CROW framework and interpretation of results.

### Setup instructions

Before you run this example, you need to install [CROW framework](https://github.com/acopar/crow) and start the container.

```sh
    git clone https://github.com/acopar/crow-example
    cd crow-example
    pip install -r requirements.txt 
```

### Use case

Before you run this example, you need to install [CROW framework](https://github.com/acopar/crow) and start the container. You can quickly set up the framework with the following commands (for details, please refer to the [CROW documentation](https://crow.readthedocs.io/).

```sh
   git clone https://github.com/acopar/crow
   cd crow
   make install
   crow-start -d
```

When the framework is running, you can run the experiment:

```sh
    python main.py
```

If you had previously installed framework to a different location, you can pass it as a variable:

```sh
   CROW_HOME='/path/to/crow' python main.py
```
