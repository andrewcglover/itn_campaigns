hipercow_finished_010425 <- read.csv("hipercow_finished_010425.csv")

hipercow_finished_010425$ISO2 <- rep(NA, dim(hipercow_finished_010425)[1])
hipercow_finished_010425$net_strategy <- rep(NA, dim(hipercow_finished_010425)[1])

for (i in 1:dim(hipercow_finished_010425)[1]) {
  test <- extract_hipercow_net_runs(hipercow_finished_010425$Job.Name[i])
  hipercow_finished_010425$ISO2[i] <- unique(test$ISO2)
  hipercow_finished_010425$net_strategy[i] <- unique(test$net_strategy)
}

#SN

test1 <- extract_hipercow_net_runs(
  "d797d5780ee6e39fa3fc9b24d2093a7c"
)
test2 <- extract_hipercow_net_runs(
  "a1daa1d054859fca0d96cd644ce478d3"
)
test3 <- extract_hipercow_net_runs(
  "aed9a5ca9fd6075e6dccf6992f63f95c"
)
test4 <- extract_hipercow_net_runs(
  "649e127e2943bb3716f37bd92f313ad6"
)
test5 <- extract_hipercow_net_runs(
  "cc8937c79a81082868473086575ecca2"
)
test6 <- extract_hipercow_net_runs(
  "c93138d5777598b8fdeb33a393b09d0f"
)
test7 <- extract_hipercow_net_runs(
  "f1bd025691b86576c5496f4f80b50414"
)
test8 <- extract_hipercow_net_runs(
  "6dcb1de996d2418776ba474cc7b3cca8"
)
test9 <- extract_hipercow_net_runs(
  "35d115896c8723af7f394b93c8cc5a4b"
)
test10 <- extract_hipercow_net_runs(
  "7e950dab27253686746638792a875838"
)

cat(
  paste(
    unique(test1$ISO2), unique(test1$net_strategy) , "\n",
    unique(test2$ISO2), unique(test2$net_strategy) , "\n",
    unique(test3$ISO2), unique(test3$net_strategy) , "\n",
    unique(test4$ISO2), unique(test4$net_strategy) , "\n",
    unique(test5$ISO2), unique(test5$net_strategy) , "\n",
    unique(test6$ISO2), unique(test6$net_strategy) , "\n",
    unique(test7$ISO2), unique(test7$net_strategy) , "\n",
    unique(test8$ISO2), unique(test8$net_strategy) , "\n",
    unique(test9$ISO2), unique(test9$net_strategy) , "\n",
    unique(test10$ISO2), unique(test10$net_strategy)
  )
)



#ML

test0 <- extract_hipercow_net_runs(
  "ff28596b4cd8286d0ecd91ed24abca27"
)
test1 <- extract_hipercow_net_runs(
  "f51ec8d06ec60c7456ed805bd7f6cebe"
)
test2 <- extract_hipercow_net_runs(
  "16466979f4301623cd6d32ad7f6fc662"
)
test3 <- extract_hipercow_net_runs(
  "1866519d8c503d8358743828f656e891"
)
test4 <- extract_hipercow_net_runs(
  "d8974fbb971c49ae7aefc84a4bbc946a"
)
test5 <- extract_hipercow_net_runs(
  "c831325fd4c3b6fca9a8fa52cad25eca"
)
test6 <- extract_hipercow_net_runs(
  "81a4bf52c6aaaa46058be1f800b2abf1"
)
test7 <- extract_hipercow_net_runs(
  "2a814ecb9117dffe82f012302acac0a8"
)
test8 <- extract_hipercow_net_runs(
  "e7dc9e7b0e34f3ef2040588303b66d18"
)
test9 <- extract_hipercow_net_runs(
  "a7da7de1dad9487c0cd7ef829fe1e851"
)
test10 <- extract_hipercow_net_runs(
  "49d87a04243f4ea3ae5d0a419feec8c0"
)

cat(
  paste(
    unique(test0$ISO2), unique(test0$net_strategy) , "\n",
    unique(test1$ISO2), unique(test1$net_strategy) , "\n",
    unique(test2$ISO2), unique(test2$net_strategy) , "\n",
    unique(test3$ISO2), unique(test3$net_strategy) , "\n",
    unique(test4$ISO2), unique(test4$net_strategy) , "\n",
    unique(test5$ISO2), unique(test5$net_strategy) , "\n",
    unique(test6$ISO2), unique(test6$net_strategy) , "\n",
    unique(test7$ISO2), unique(test7$net_strategy) , "\n",
    unique(test8$ISO2), unique(test8$net_strategy) , "\n",
    unique(test9$ISO2), unique(test9$net_strategy) , "\n",
    unique(test10$ISO2), unique(test10$net_strategy)
  )
)







test1 <- extract_hipercow_net_runs(
  "f51ec8d06ec60c7456ed805bd7f6cebe"
)
test2 <- extract_hipercow_net_runs(
  "16466979f4301623cd6d32ad7f6fc662"
)
test3 <- extract_hipercow_net_runs(
  "1866519d8c503d8358743828f656e891"
)
test4 <- extract_hipercow_net_runs(
  "ff28596b4cd8286d0ecd91ed24abca27"
)
test5 <- extract_hipercow_net_runs(
  "93563db696992a33c040e515db0cd386"
)
test6 <- extract_hipercow_net_runs(
  "d10571627cfd4347ac03f8703b40ca85"
)
test7 <- extract_hipercow_net_runs(
  "c037008ea3ed8309c8df1532a80c3e5b"
)
test8 <- extract_hipercow_net_runs(
  "2ed59b4080a09fd4b94ba878fd737620"
)
test9 <- extract_hipercow_net_runs(
  "c24a9303bacbabe97dad8572117d37ce"
)
test10 <- extract_hipercow_net_runs(
  "49d87a04243f4ea3ae5d0a419feec8c0"
)

cat(
  paste(
    unique(test1$ISO2), unique(test1$net_strategy) , "\n",
    unique(test2$ISO2), unique(test2$net_strategy) , "\n",
    unique(test3$ISO2), unique(test3$net_strategy) , "\n",
    unique(test4$ISO2), unique(test4$net_strategy) , "\n",
    unique(test5$ISO2), unique(test5$net_strategy) , "\n",
    unique(test6$ISO2), unique(test6$net_strategy) , "\n",
    unique(test7$ISO2), unique(test7$net_strategy) , "\n",
    unique(test8$ISO2), unique(test8$net_strategy) , "\n",
    unique(test9$ISO2), unique(test9$net_strategy) , "\n",
    unique(test10$ISO2), unique(test10$net_strategy)
  )
)