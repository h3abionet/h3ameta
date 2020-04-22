

params.message = "hello"

data_ch = Channel.from(params.message)

process banner {
  input:
    val(message) from data_ch
  output:
    file(result) into banner_ch
  script:
    result = "banner.txt"
    """
    banner ${params.message} >  $result
    """
}


process countRows {
  input:
    file(result) from banner_ch
  output:
    stdout col_ch
  script:
  """
    cut -b 10 $result > col.txt
    grep -c \\#  col.txt
  """
}

col_ch.subscribe { println it }
