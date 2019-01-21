

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
  validExitStatus 0, 1
  input:
    file(result) from banner_ch
  each col from 1..80
  output:
    stdout col_ch
  script:
  """
    cut -b $col $result > col.txt
    grep -c \\#  col.txt
  """
}

process combineResults {
  input:
  val(res) from col_ch.toList()
  output:
    stdout comb_ch
  script:
    data = res*.trim().join(" ")
    """
     echo  $data
   """
}
comb_ch.subscribe { println it }
