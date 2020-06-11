import click


@click.group()
def main():
	pass


@main.group()
def test1():
	pass


@main.group()
def test2():
	pass


@test1.command()
@click.option('--sub1')
def sub1():
	pass


@test2.command()
@click.option('--sub2')
def sub2():
	pass


if __name__ == '__main__':
	main()