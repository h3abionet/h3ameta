

params.message = "hello"

data_ch = Channel.from(params.message)


data_ch.subscribe { println it }
